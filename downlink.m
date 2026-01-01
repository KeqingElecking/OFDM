% run_nbiot_downlink.m
clear; clc; close all;

%% 1. Parameters
ENABLE_EQ = true;
N_fft = 64; N_cp = 16; 
N_sym = 14*2; % 2 Subframes
CellID = 10; M = 4;
SNR_dB = 20;

%% 2. Downlink Tx
disp('--- Generating Downlink Signal (OFDMA) ---');

% 2a. Dynamic Pilot Calculation (Standard NRS)
pilot_symbols = [];
for s = 1:N_sym
    if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
        pilot_symbols = [pilot_symbols, s];
    end
end

% 2b. Generate Data
% Calculate exact capacity
num_pilots = 0;
for s = 1:N_sym
    if ismember(s, pilot_symbols), num_pilots = num_pilots + 2; end
end
n_bits = (12*N_sym - num_pilots) * 2;
tx_bits = randi([0 1], n_bits, 1);

% 2c. Map & Modulate
% Note: Ensure your 'nbiot_std_mapper' uses the dynamic pilot logic too!
% Or define specific mapper inside script if not modifying file.
[grid_dl, p_mask] = nbiot_std_mapper(tx_bits, N_sym, CellID); 
tx_wave = nbiot_ofdm_modulator(grid_dl, N_fft, N_cp);

%% 3. Channel
[rx_wave, H_true] = nbiot_mcm_channel(tx_wave, N_fft, N_cp, SNR_dB);

%% 4. Downlink Rx
[~, rx_grid_raw] = nbiot_rx_simple(rx_wave, N_fft, N_cp, N_sym, CellID);

if ENABLE_EQ
    [rx_grid_final, H_est] = nbiot_equalizer(rx_grid_raw, CellID);
else
    rx_grid_final = rx_grid_raw;
end

%% 5. Demodulate & BER
rx_data = rx_grid_final(~p_mask);
rx_bits = qamdemod(rx_data, M, 'OutputType', 'bit', 'UnitAveragePower', true);

len = min(length(tx_bits), length(rx_bits));
[~, ber] = biterr(tx_bits(1:len), rx_bits(1:len));

fprintf('Downlink BER: %.4f\n', ber);

figure('Name', 'Downlink Analysis');
subplot(1,2,1);
plot(rx_data, 'r.'); title('DL Constellation'); axis square; grid on;
subplot(1,2,2);
imagesc(abs(rx_grid_final)); title('Equalized Grid');