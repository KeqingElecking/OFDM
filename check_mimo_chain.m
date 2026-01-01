% CHECK_MIMO_CHAIN.m
% Verification for NB-IoT 2x2 MIMO Transmit Diversity
clear; clc; close all;

%% 1. Parameters
N_fft = 64; N_cp = 16; N_sym = 14; 
CellID = 5; % Try changing to see pilot shifts
M = 4;

%% 2. Generate Data
% Need enough bits for 2 antennas? No, SFBC sends SAME bits redundancy.
% So we generate bits for roughly 1 stream capacity.
[tx_bits, ~] = gen_ofdm_data(12, N_sym, M, 1); 

%% 3. MIMO Mapping (SFBC)
[grid_ant0, grid_ant1] = nbiot_mimo_mapper(tx_bits, N_sym, CellID);

% VISUALIZATION: Compare Antennas
figure('Name', 'MIMO Grids');
subplot(1,2,1);
imagesc(abs(grid_ant0)); title('Antenna 0 (Pilots + Data)');
xlabel('Sym'); ylabel('Subcarrier'); colorbar;

subplot(1,2,2);
imagesc(abs(grid_ant1)); title('Antenna 1 (Shifted Pilots + SFBC Data)');
xlabel('Sym'); ylabel('Subcarrier'); colorbar;

%% 4. OFDM Modulation (Per Antenna)
% Use the SAME modulator function for both streams
tx_wave_ant0 = nbiot_ofdm_modulator(grid_ant0, N_fft, N_cp);
tx_wave_ant1 = nbiot_ofdm_modulator(grid_ant1, N_fft, N_cp);

disp('MIMO Generation Complete.');
disp(['Antenna 0 Power: ' num2str(mean(abs(tx_wave_ant0).^2))]);
disp(['Antenna 1 Power: ' num2str(mean(abs(tx_wave_ant1).^2))]);

%% 5. What's Next? (The Channel)
% To simulate the reception, you would now mix these:
% Rx = (H0 * tx_ant0) + (H1 * tx_ant1) + Noise