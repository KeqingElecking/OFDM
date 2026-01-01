% run_nbiot_uplink.m
clear; clc; close all;

%% 1. Parameters
N_fft = 64; N_cp = 16; 
N_sym = 28; % 2 Subframes (4 slots)
M = 4; % QPSK
SNR_dB = 25;

%% 2. Uplink Tx
disp('--- Generating Uplink Signal (SC-FDMA) ---');
% Generate enough bits (approx)
[tx_bits, ~] = gen_ofdm_data(12, N_sym, M, 1);
[grid_ul, p_mask] = nbiot_uplink_mapper(tx_bits, N_sym, M);

% Modulate (DFT -> Map -> IFFT)
tx_wave = nbiot_scfdma_modulator(grid_ul, N_fft, N_cp);

%% 3. Channel
disp('--- Applying Channel ---');
[rx_wave, H_true] = nbiot_mcm_channel(tx_wave, N_fft, N_cp, SNR_dB);

%% 4. Uplink Rx
disp('--- Receiving Uplink Signal ---');
[rx_bits, rx_grid_eq] = nbiot_uplink_rx(rx_wave, N_fft, N_cp, N_sym, M);

%% 5. Analysis
% We need to align Tx bits with Rx bits (handling truncation in mapper)
len = min(length(tx_bits), length(rx_bits));
[~, ber] = biterr(tx_bits(1:len), rx_bits(1:len));

fprintf('Uplink BER: %.4f\n', ber);

% Visualization
figure('Name', 'Uplink Analysis');
subplot(1,2,1);
% Plot data symbols only
data_syms = rx_grid_eq(~p_mask);
plot(data_syms, 'b.'); title(['UL Constellation (SNR=' num2str(SNR_dB) ')']);
axis square; grid on;

subplot(1,2,2);
% Visualize the Equalized Grid (Magnitude)
imagesc(abs(rx_grid_eq)); colorbar;
title('Equalized Rx Grid (Time Domain Symbols)');
xlabel('Symbol'); ylabel('Subcarrier');