% run_all_mimo_plots_parallel.m
% Parallelized NB-IoT MIMO Simulation
% Requires: Parallel Computing Toolbox
clear; clc; close all;

%% 1. Configuration
N_fft = 64; N_cp = 32; 
N_sym = 10000;      % High symbol count is fine with parallel processing
CellID = 2; M = 4;
SNR_vec = 0:1:25;   % Loop for BER Curve

% Pre-allocate Output Arrays
BER_DL_SFBC = zeros(size(SNR_vec));
BER_UL_MRC  = zeros(size(SNR_vec));

% Prepare a container for plotting data (to handle parfor limitations)
plot_data = cell(1, length(SNR_vec));

% Setup Parallel Pool (Start workers if not already running)
if isempty(gcp('nocreate'))
    parpool; % Starts the default parallel pool
end
disp(['Running Simulation on ' num2str(length(SNR_vec)) ' SNR points in parallel...']);

%% 2. Parallel Simulation Loop
% We use 'parfor' here. Note that 'i' is the loop variable.
parfor i = 1:length(SNR_vec)
    snr = SNR_vec(i);
    % fprintf inside parfor might appear out of order, but helps track progress
    fprintf('Worker processing SNR: %d dB...\n', snr);
    
    % ==========================================================
    % --- DOWNLINK (2 Tx -> 1 Rx, SFBC) ---
    % ==========================================================
    % Note: Functions like gen_ofdm_data must be on the path
    [tx_bits_dl, ~] = gen_ofdm_data(12, N_sym, M, 1);
    [grid_ant0, grid_ant1] = nbiot_mimo_mapper(tx_bits_dl, N_sym, CellID);
    
    tx_w0 = nbiot_ofdm_modulator(grid_ant0, N_fft, N_cp);
    tx_w1 = nbiot_ofdm_modulator(grid_ant1, N_fft, N_cp);
    
    % Run Channel
    [rx_dl, ~, H_struct_dl] = nbiot_mimo_channel(tx_w0, tx_w1, N_fft, N_cp, snr);
    
    % Receiver
    [rx_bits_dl, grid_eq_dl] = nbiot_dl_mimo_rx(rx_dl, N_fft, N_cp, N_sym, CellID);
    
    % Error Count
    len_dl = min(length(tx_bits_dl), length(rx_bits_dl));
    [~, ber_dl] = biterr(tx_bits_dl(1:len_dl), rx_bits_dl(1:len_dl));
    BER_DL_SFBC(i) = ber_dl; % Sliced output assignment
    
    % ==========================================================
    % --- UPLINK (1 Tx -> 2 Rx, MRC) ---
    % ==========================================================
    [tx_bits_ul, ~] = gen_ofdm_data(12, N_sym, M, 1);
    [grid_ul, ~] = nbiot_uplink_mapper(tx_bits_ul, N_sym, M);
    tx_w_ul = nbiot_scfdma_modulator(grid_ul, N_fft, N_cp);
    
    % 1x2 Channel (Pass zeros for 2nd Tx antenna)
    [rx_ul_1, rx_ul_2, ~] = nbiot_mimo_channel(tx_w_ul, zeros(size(tx_w_ul)), N_fft, N_cp, snr);
    
    % Receiver
    [rx_bits_ul, grid_eq_ul] = nbiot_ul_mimo_rx(rx_ul_1, rx_ul_2, N_fft, N_cp, N_sym, M);
    
    % Error Count
    len_ul = min(length(tx_bits_ul), length(rx_bits_ul));
    [~, ber_ul] = biterr(tx_bits_ul(1:len_ul), rx_bits_ul(1:len_ul));
    BER_UL_MRC(i) = ber_ul; % Sliced output assignment
    
    % ==========================================================
    % --- SAVE DATA FOR PLOTTING (Conditional) ---
    % ==========================================================
    % We cannot save to a single variable in parfor. We save to the cell array.
    if snr == 15
        s = struct();
        s.grid_dl = grid_eq_dl;
        s.grid_ul = grid_eq_ul;
        s.H_dl_path1 = H_struct_dl.H11;
        plot_data{i} = s;
    end
end

disp('Simulation Complete. Extracting Plot Data...');

%% 3. Data Extraction (Post-Processing)
% Find the index where SNR was 15 and extract the struct
idx_plot = find(SNR_vec == 15, 1);

if ~isempty(idx_plot) && ~isempty(plot_data{idx_plot})
    saved_grid_dl = plot_data{idx_plot}.grid_dl;
    saved_grid_ul = plot_data{idx_plot}.grid_ul;
    saved_H_dl    = plot_data{idx_plot}.H_dl_path1;
else
    warning('SNR=15dB data not found. Visualization will be skipped or empty.');
    saved_grid_dl = []; saved_grid_ul = []; saved_H_dl = [];
end

%% 4. Visualization (Standard Plotting Code)

% Figure 1: BER Performance
figure('Name', 'NB-IoT MIMO Performance', 'Position', [100 100 600 400]);
semilogy(SNR_vec, BER_DL_SFBC, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(SNR_vec, BER_UL_MRC, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title(['NB-IoT MIMO (Parallel Run, N_{sym}=' num2str(N_sym) ')']);
legend('Downlink (SFBC 2x1)', 'Uplink (MRC 1x2)');
ylim([1e-6 1]); % Adjusted lower limit for high N_sym

% Check if we have plotting data before plotting figs 2 & 3
if ~isempty(saved_grid_dl)
    % Figure 2: Constellations at SNR = 15dB
    figure('Name', 'Constellation Diagrams (15dB)', 'Position', [750 100 800 400]);
    subplot(1,2,1);
    plot(saved_grid_dl(:), 'b.'); hold on;
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
    [L_taps, n_syms_plotting] = size(saved_H_dl);
    % Limit plotting to first 100 symbols to avoid massive heavy plot
    n_syms_plotting = min(n_syms_plotting, 200); 
    
    H_freq_response = zeros(N_fft, n_syms_plotting);
    for s = 1:n_syms_plotting
        taps = saved_H_dl(:, s);
        H_freq_response(:, s) = fft(taps, N_fft); 
    end
    H_freq_centered = fftshift(H_freq_response, 1);
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
end