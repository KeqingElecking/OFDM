% CHECK_EQUALIZER_TOGGLE_DYNAMIC.m
clear; clc;

%% 1. Parameters
% --- CONFIGURATION ---
ENABLE_EQUALIZER = true;  % Set to false to disable Equalization
% ---------------------
N_fft = 64; N_cp = 32; 
N_sym = 28;  % CHANGED: Try 28 (2 Subframes) to test the scaling!
CellID = 2; M = 4; k = 2;
SNR_dB = 25; % Increased SNR slightly for clearer constellation 

%% 2. Data & Tx
% --- UPDATED: Dynamic Pilot Calculation ---
% Standard NB-IoT/LTE: Pilots are in symbol indices 6, 7 of every slot.
% (Slot = 7 symbols). So indices are 6, 7, 13, 14, 20, 21...
pilot_syms = [];
for s = 1:N_sym
    % Check if 's' is the 6th or 7th symbol in a slot (modulo 7)
    % mod(s-1, 7) gives 0..6. We want indices 5 and 6 (which map to 6th and 7th)
    if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
        pilot_syms = [pilot_syms, s];
    end
end
% ------------------------------------------

% Calculate Capacity
num_pilots = 0;
for s=1:N_sym
    if ismember(s,pilot_syms)
        % Standard NRS has 2 pilots per symbol (Ports 2000/2001 logic)
        num_pilots = num_pilots + 2; 
    end
end
data_bits_len = (12*N_sym - num_pilots) * k;

tx_bits = randi([0 1], data_bits_len, 1);

% Transmit Chain
% NOTE: Ensure your nbiot_std_mapper handles N_sym > 14 correctly!
[grid_tx, ~] = nbiot_std_mapper(tx_bits, N_sym, CellID);
tx_wave = nbiot_ofdm_modulator(grid_tx, N_fft, N_cp);

%% 3. Channel
disp('Applying Channel...');
[rx_wave, H_true] = nbiot_mcm_channel(tx_wave, N_fft, N_cp, SNR_dB);

%% 4. Receiver Front-End
[~, rx_grid_raw] = nbiot_rx_simple(rx_wave, N_fft, N_cp, N_sym, CellID);

%% 5. Equalization
if ENABLE_EQUALIZER
    disp('Equalizer Status: ON');
    [rx_grid_final, H_est] = nbiot_equalizer(rx_grid_raw, CellID);
else
    disp('Equalizer Status: OFF (Bypassing)');
    rx_grid_final = rx_grid_raw;
    H_est = ones(size(rx_grid_raw));
end

%% 6. Demodulation
% Create Mask (Dynamic based on updated pilot_syms)
mask = false(12, N_sym);
for s=1:N_sym
    if ismember(s, pilot_syms)
        p = [1 7] + mod(CellID,6); 
        p(p>12) = p(p>12)-12;
        mask(p,s) = true;
    end
end

% A. Demodulate RAW
raw_syms = rx_grid_raw(~mask);
bits_raw = qamdemod(raw_syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);
% Handle truncation if sizes mismatch slightly due to mapping logic
len = min(length(tx_bits), length(bits_raw));
[~, ber_raw] = biterr(tx_bits(1:len), bits_raw(1:len));

% B. Demodulate FINAL
final_syms = rx_grid_final(~mask);
bits_final = qamdemod(final_syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);
[~, ber_final] = biterr(tx_bits(1:len), bits_final(1:len));

%% 7. Visualization
figure('Position', [100 100 1000 400]);

% Plot 1: Channel Estimation Accuracy
subplot(1,3,1);
plot(abs(H_true(1, :)), 'g-', 'LineWidth', 2); hold on;
plot(abs(H_est(1, :)), 'r--', 'LineWidth', 2);
title('Channel Est (Subcarrier 1)');
xlabel('Symbol'); 
if ENABLE_EQUALIZER, legend('True Channel', 'Estimated H');
else, legend('True Channel', 'Assume Flat (Eq Off)'); end
grid on;

% Plot 2: Constellations
subplot(1,3,2);
plot(raw_syms, 'r.'); hold on;
plot(final_syms, 'b.'); 
title('Constellation');
legend(['Raw (BER ' num2str(ber_raw, '%.2f') ')'], ...
       ['Final (BER ' num2str(ber_final, '%.4f') ')']);
axis square; grid on;

% Plot 3: 3D View
subplot(1,3,3);
surf(abs(H_est)); shading interp;
title('Estimated Channel Magnitude');
xlabel('Symbol'); ylabel('Subcarrier'); zlabel('|H|');
view(-45, 30);
if ~ENABLE_EQUALIZER, zlim([0 2]); end