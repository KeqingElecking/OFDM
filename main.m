% run_all_mimo_plots.m
% Comprehensive NB-IoT MIMO Simulation with Visualization
clear; clc; close all;

%% 1. Configuration
N_fft = 64; N_cp = 32; 
N_sym = 10000; % Increased symbols for smoother BER curves
CellID = 2; M = 4;
SNR_vec = 0:1:25; % Loop for BER Curve

% Pre-allocate BER arrays
BER_DL_SFBC = zeros(size(SNR_vec));
BER_UL_MRC  = zeros(size(SNR_vec));

disp('Running Simulation Loop (this may take a moment)...');

%% 2. Simulation Loop (BER vs SNR)
for i = 1:length(SNR_vec)
    snr = SNR_vec(i);
    fprintf('  Simulating SNR: %d dB...\n', snr);
    
    % --- DOWNLINK (2 Tx -> 1 Rx, SFBC) ---
    [tx_bits_dl, ~] = gen_ofdm_data(12, N_sym, M, 1);
    [grid_ant0, grid_ant1] = nbiot_mimo_mapper(tx_bits_dl, N_sym, CellID);
    
    tx_w0 = nbiot_ofdm_modulator(grid_ant0, N_fft, N_cp);
    tx_w1 = nbiot_ofdm_modulator(grid_ant1, N_fft, N_cp);
    
    [rx_dl, ~, H_struct_dl] = nbiot_mimo_channel(tx_w0, tx_w1, N_fft, N_cp, snr);
    
    [rx_bits_dl, grid_eq_dl] = nbiot_dl_mimo_rx(rx_dl, N_fft, N_cp, N_sym, CellID);
    
    len_dl = min(length(tx_bits_dl), length(rx_bits_dl));
    [~, BER_DL_SFBC(i)] = biterr(tx_bits_dl(1:len_dl), rx_bits_dl(1:len_dl));
    
    % --- UPLINK (1 Tx -> 2 Rx, MRC) ---
    [tx_bits_ul, ~] = gen_ofdm_data(12, N_sym, M, 1);
    [grid_ul, ~] = nbiot_uplink_mapper(tx_bits_ul, N_sym, M);
    tx_w_ul = nbiot_scfdma_modulator(grid_ul, N_fft, N_cp);
    
    % 1x2 Channel (Pass zeros for 2nd Tx antenna)
    [rx_ul_1, rx_ul_2, ~] = nbiot_mimo_channel(tx_w_ul, zeros(size(tx_w_ul)), N_fft, N_cp, snr);
    
    [rx_bits_ul, grid_eq_ul] = nbiot_ul_mimo_rx(rx_ul_1, rx_ul_2, N_fft, N_cp, N_sym, M);
    
    len_ul = min(length(tx_bits_ul), length(rx_bits_ul));
    [~, BER_UL_MRC(i)] = biterr(tx_bits_ul(1:len_ul), rx_bits_ul(1:len_ul));
    
    % Save data for plotting constellations at specific SNR (e.g., 15dB)
    if snr == 15
        saved_grid_dl = grid_eq_dl;
        saved_grid_ul = grid_eq_ul;
        saved_H_dl = H_struct_dl.H11; % Just grab one path for visualization
    end
end

%% 3. Visualization

% Figure 1: BER Performance
figure('Name', 'NB-IoT MIMO Performance', 'Position', [100 100 600 400]);
semilogy(SNR_vec, BER_DL_SFBC, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(SNR_vec, BER_UL_MRC, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title('NB-IoT MIMO Performance Comparison');
legend('Downlink (SFBC 2x1)', 'Uplink (MRC 1x2)');
ylim([1e-5 1]);

% Figure 2: Constellations at SNR = 15dB
figure('Name', 'Constellation Diagrams (15dB)', 'Position', [750 100 800 400]);

subplot(1,2,1);
plot(saved_grid_dl(:), 'b.'); hold on;
% Ideal QPSK points
qpsk_pts = [1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2);
plot(qpsk_pts, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
title('Downlink SFBC Equalized Symbols');
axis([-2 2 -2 2]); axis square; grid on;
legend('Rx Symbols', 'Ideal QPSK');

subplot(1,2,2);
plot(saved_grid_ul(:), 'm.'); hold on;
plot(qpsk_pts, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
title('Uplink MRC Equalized Symbols');
axis([-2 2 -2 2]); axis square; grid on;
legend('Rx Symbols', 'Ideal QPSK');

% Figure 3: Channel Physics (One Path)
% We reshape the 1D channel vector back to the grid to see time/freq fading
% The channel history 'saved_H_dl' contains the Impulse Response (Taps)
% Dimensions: [L_taps x N_sym] (e.g., 6 x 100)
% We must apply FFT to convert Taps -> Frequency Response

[L_taps, n_syms_plotting] = size(saved_H_dl);
H_freq_response = zeros(N_fft, n_syms_plotting);

for s = 1:n_syms_plotting
    taps = saved_H_dl(:, s);
    % FFT of taps gives the frequency response for that OFDM symbol
    % We FFT to N_fft size (64) to see the full band
    H_freq_response(:, s) = fft(taps, N_fft); 
end

% Shift zero frequency to center to match our subcarrier mapping
H_freq_centered = fftshift(H_freq_response, 1);

% Extract only the 12 active subcarriers where our data lives
% In Modulator: zeros_left = floor((64-12)/2) = 26. Start index = 27.
center_idx = floor((N_fft - 12)/2) + 1; 
active_indices = center_idx : (center_idx + 11);

H_active = H_freq_centered(active_indices, :);

figure('Name', 'Channel Fading Profile', 'Position', [400 550 600 400]);
surf(1:n_syms_plotting, 1:12, abs(H_active)); 
shading interp;
title('Fading Channel Magnitude (Active Subcarriers)');
xlabel('OFDM Symbol (Time)'); 
ylabel('Subcarrier Index (Freq)');
zlabel('Magnitude |H|');
colorbar;
view(-45, 30);