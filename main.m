% run_all_mimo_plots.m
% Comprehensive NB-IoT MIMO Simulation with Visualization
% Features: Frozen Channel, MMSE 2D, Dual Channel Physics Plot
clear; clc; close all;

%% 1. Configuration
N_fft = 64; N_cp = 32; 
N_sym = 10000; 
CellID = 2; M = 4;
SNR_vec = 0:2:24; 

% Pre-allocate BER arrays
BER_DL_SFBC = zeros(size(SNR_vec));
BER_UL_MRC  = zeros(size(SNR_vec));

% Frozen Channel Seed
MASTER_SEED = 12345; 

disp('Running Simulation Loop (this may take a moment)...');

%% 2. Simulation Loop
for i = 1:length(SNR_vec)
    snr = SNR_vec(i);
    fprintf('  Simulating SNR: %d dB...\n', snr);
    
    % Force identical channel physics for every SNR point
    rng(MASTER_SEED); 
    
    % --- DOWNLINK (2 Tx -> 1 Rx, SFBC) ---
    [tx_bits_dl, ~] = gen_ofdm_data(12, N_sym, M, 1);
    [grid_ant0, grid_ant1] = nbiot_mimo_mapper(tx_bits_dl, N_sym, CellID);
    
    tx_w0 = nbiot_ofdm_modulator(grid_ant0, N_fft, N_cp);
    tx_w1 = nbiot_ofdm_modulator(grid_ant1, N_fft, N_cp);
    
    [rx_dl, ~, H_struct_dl] = nbiot_mimo_channel(tx_w0, tx_w1, N_fft, N_cp, snr);
    
    [rx_bits_dl, grid_eq_dl] = nbiot_dl_mimo_rx(rx_dl, N_fft, N_cp, N_sym, CellID, snr);
    
    len_dl = min(length(tx_bits_dl), length(rx_bits_dl));
    [~, BER_DL_SFBC(i)] = biterr(tx_bits_dl(1:len_dl), rx_bits_dl(1:len_dl));
    
    % --- UPLINK (1 Tx -> 2 Rx, MRC) ---
    [tx_bits_ul, ~] = gen_ofdm_data(12, N_sym, M, 1);
    [grid_ul, ~] = nbiot_uplink_mapper(tx_bits_ul, N_sym, M);
    tx_w_ul = nbiot_scfdma_modulator(grid_ul, N_fft, N_cp);
    
    % Capture Uplink Channel Struct (H_struct_ul)
    [rx_ul_1, rx_ul_2, H_struct_ul] = nbiot_mimo_channel(tx_w_ul, zeros(size(tx_w_ul)), N_fft, N_cp, snr);
    
    [rx_bits_ul, grid_eq_ul] = nbiot_ul_mimo_rx(rx_ul_1, rx_ul_2, N_fft, N_cp, N_sym, M, snr);
    
    len_ul = min(length(tx_bits_ul), length(rx_bits_ul));
    [~, BER_UL_MRC(i)] = biterr(tx_bits_ul(1:len_ul), rx_bits_ul(1:len_ul));
    
    % Save data for visualization at 14dB
    if snr == 14
        saved_grid_dl = grid_eq_dl;
        saved_grid_ul = grid_eq_ul;
        
        % Save Channel Impulse Responses
        saved_H_dl = H_struct_dl.H11; % Downlink Path (Tx1 -> Rx1)
        saved_H_ul = H_struct_ul.H11; % Uplink Path (Tx1 -> Rx1)
    end
end

%% 3. Visualization

% Figure 1: BER
figure('Name', 'NB-IoT MIMO Performance', 'Position', [100 100 600 400]);
semilogy(SNR_vec, BER_DL_SFBC, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(SNR_vec, BER_UL_MRC, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('NB-IoT MIMO Performance Comparison');
legend('Downlink (SFBC 2x1)', 'Uplink (MRC 1x2)'); ylim([1e-5 1]);

% Figure 2: Constellations (14dB)
figure('Name', 'Constellations (14dB)', 'Position', [750 100 800 400]);
subplot(1,2,1);
plot(saved_grid_dl(:), 'b.'); hold on;
plot([1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2), 'ro', 'LineWidth', 2);
title('Downlink SFBC Equalized'); axis square; grid on;

subplot(1,2,2);
plot(saved_grid_ul(:), 'm.'); hold on;
plot([1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2), 'ro', 'LineWidth', 2);
title('Uplink MRC Equalized'); axis square; grid on;

% --- FIGURE 3: DUAL CHANNEL PHYSICS (DL & UL) ---
figure('Name', 'Channel Fading Profile', 'Position', [400 550 600 700]);

% Prepare Axes
[~, n_syms_plot] = size(saved_H_dl);
subcarrier_spacing = 15e3; 
freq_axis_khz = ((-6:5) * subcarrier_spacing) / 1000; 
center_idx = floor((N_fft - 12)/2) + 1; 
active_idx = center_idx : (center_idx + 11);

% --- TOP PLOT: DOWNLINK CHANNEL ---
H_freq_dl = zeros(N_fft, n_syms_plot);
for s = 1:n_syms_plot, H_freq_dl(:, s) = fft(saved_H_dl(:, s), N_fft); end
H_centered_dl = fftshift(H_freq_dl, 1);
H_active_dl = H_centered_dl(active_idx, :);

subplot(2,1,1);
surf(1:n_syms_plot, freq_axis_khz, abs(H_active_dl)); 
shading interp; colorbar; view(-45, 30);
title('Downlink Channel Magnitude (Tx1 \rightarrow Rx1)');
ylabel('Freq Offset (kHz)'); zlabel('|H|');
axis tight;

% --- BOTTOM PLOT: UPLINK CHANNEL ---
H_freq_ul = zeros(N_fft, n_syms_plot);
for s = 1:n_syms_plot, H_freq_ul(:, s) = fft(saved_H_ul(:, s), N_fft); end
H_centered_ul = fftshift(H_freq_ul, 1);
H_active_ul = H_centered_ul(active_idx, :);

subplot(2,1,2);
surf(1:n_syms_plot, freq_axis_khz, abs(H_active_ul)); 
shading interp; colorbar; view(-45, 30);
title('Uplink Channel Magnitude (Tx1 \rightarrow Rx1)');
xlabel('OFDM Symbol (Time)'); ylabel('Freq Offset (kHz)'); zlabel('|H|');
axis tight;