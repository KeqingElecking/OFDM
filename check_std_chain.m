% CHECK_STD_CHAIN.m
% Verification for Standard NB-IoT Downlink Chain (NRS Pilots)
clear; clc; close all;

%% 1. Parameters
N_fft = 64;
N_cp = 16;
N_sym = 14;      % 1 Subframe
CellID = 2;      % Try changing this (0-5) to see pilots move!
M = 4;           % QPSK

%% 2. Generate Data
% Using your uploaded function
[tx_bits, ~] = gen_ofdm_data(12, N_sym, M, 1); 

%% 3. Map to Grid (Conventional NRS)
[grid_std, pilot_mask] = nbiot_std_mapper(tx_bits, N_sym, CellID);

% VISUALIZATION 1: Resource Grid
figure('Name', 'Standard NB-IoT Grid');
imagesc(abs(grid_std));
colorbar;
title(['Standard NRS Pattern (Cell ID: ' num2str(CellID) ')']);
xlabel('OFDM Symbols'); ylabel('Subcarriers');
set(gca, 'YTick', 1:12);
set(gca, 'XTick', 1:14);

% Overlay Pilots
hold on;
[r, c] = find(pilot_mask);
plot(c, r, 'wo', 'LineWidth', 2, 'MarkerSize', 8);
legend('NRS Pilots');
grid on;

%% 4. OFDM Modulation
% Using your uploaded function
tx_waveform = nbiot_ofdm_modulator(grid_std, N_fft, N_cp);

%% 5. Output Analysis
figure('Name', 'Tx Signal Analysis');

% Time Domain
subplot(2,1,1);
plot(real(tx_waveform));
xline(N_fft+N_cp, 'r--', 'Symbol 1 End');
title('Time Domain Signal (Real)');
grid on;

% Frequency Domain
subplot(2,1,2);
% Take FFT of the whole signal to see occupied bandwidth
L = length(tx_waveform);
f = (-L/2 : L/2 - 1) * (15e3 / N_fft); % approx freq axis
spectrum = fftshift(fft(tx_waveform));
plot(abs(spectrum));
title('Full Spectrum (Check Bandwidth)');
xline([L/2 - 6*L/N_fft, L/2 + 6*L/N_fft], 'r--', 'Expected 180kHz');
grid on;