% VERIFY_MCM_CHANNEL.m
clear; clc;

%% 1. Parameters
N_fft = 64; 
N_cp = 16; 
N_sym = 20; % Simulate 20 symbols to see fading changes
M = 4;      % QPSK
snr_dB = 25;

%% 2. Generate Tx Signal (Using your previous functions)
[tx_bits, ~] = gen_ofdm_data(12, N_sym, M, 1);
[grid, ~] = nbiot_std_mapper(tx_bits, N_sym, 1);
tx_wave = nbiot_ofdm_modulator(grid, N_fft, N_cp);

%% 3. Apply New MCM Channel
% This replaces the 'conv' loop in your old code
[rx_wave, H_history] = nbiot_sos_channel(tx_wave, N_fft, N_cp, snr_dB);

%% 4. Visualization
figure;

% Plot 1: Channel Magnitude variation over time
subplot(2,1,1);
% Plot magnitude of the 1st tap (Main path) over symbols
stem(abs(H_history(1,:)), 'LineWidth', 1.5);
title('Tap 1 Magnitude Variation (Doppler Effect)');
xlabel('OFDM Symbol Index'); ylabel('Magnitude |h|');
grid on;

% Plot 2: Received Waveform vs Transmitted
subplot(2,1,2);
plot(real(tx_wave(1:200)), 'b'); hold on;
plot(real(rx_wave(1:200)), 'r--');
title('Tx vs Rx Waveform (First 200 samples)');
legend('Transmitted', 'Received (Faded+Noise)');
grid on;

disp('Channel Simulation Complete.');