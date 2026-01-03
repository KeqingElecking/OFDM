% compare_coding_metrics_final.m
% Compares Uncoded vs Hamming vs Repetition coding in NB-IoT
% - Uses "Frozen Channel" for fair comparison
% - Calculates Throughput and Overhead metrics
clear; clc; close all;

%% 1. Configuration
N_fft = 64; N_cp = 16; 
N_sym = 10000;       % Simulation Duration (Symbols)
CellID = 2; M = 4;  % QPSK
SNR_vec = 0:1:10;   % SNR Sweep

% Master Seed to freeze the physics
MASTER_SEED = 1000; 

% Storage
BER_Uncoded = zeros(size(SNR_vec));
BER_Hamming = zeros(size(SNR_vec));
BER_Repetit = zeros(size(SNR_vec));

% --- CAPACITY CALCULATION (Dynamic) ---
% We calculate how many bits fit in the grid dynamically
% (Must match the logic in your fixed nbiot_std_mapper)
pilot_count = 0;
for s = 1:N_sym
    if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
        pilot_count = pilot_count + 2; 
    end
end
grid_capacity_res = (12 * N_sym) - pilot_count;
grid_capacity_bits = grid_capacity_res * 2; 

% Throughput Timing: 14 symbols = 1ms
sim_duration_sec = N_sym * (1e-3 / 14);

fprintf('--- SETUP ---\n');
fprintf('Duration: %.4f seconds\n', sim_duration_sec);
fprintf('Capacity: %d bits\n', grid_capacity_bits);
disp('Starting Simulation Loop...');

%% 2. Simulation Loop
for i = 1:length(SNR_vec)
    snr = SNR_vec(i);
    fprintf('Simulating SNR: %d dB...\n', snr);
    
    % =======================================================
    % SCHEME A: UNCODED
    % =======================================================
    tx_data = randi([0 1], grid_capacity_bits, 1);
    
    [grid, ~] = nbiot_std_mapper(tx_data, N_sym, CellID);
    tx_wave = nbiot_ofdm_modulator(grid, N_fft, N_cp);
    
    rng(MASTER_SEED + i*100); % FREEZE CHANNEL
    rx_wave = nbiot_mcm_channel(tx_wave, N_fft, N_cp, snr);
    
    [~, rx_grid] = nbiot_rx_simple(rx_wave, N_fft, N_cp, N_sym, CellID);
    [~, H_est] = nbiot_equalizer(rx_grid, CellID);
    
    % Manual Demod (to ensure mask alignment)
    p_mask = (abs(H_est) == 0); % Or regenerate mask locally if equalizer returns 0s
    % Robust Mask Regeneration:
    p_mask = false(12, N_sym);
    for s=1:N_sym, if mod(s-1,7)==5||mod(s-1,7)==6, p_subs=[1,7]+mod(CellID,6); p_subs(p_subs>12)=p_subs(p_subs>12)-12; p_mask(p_subs,s)=true; end; end
    
    rx_val = rx_grid(~p_mask) ./ H_est(~p_mask);
    rx_bits = qamdemod(rx_val, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    
    L = min(length(tx_data), length(rx_bits));
    [~, BER_Uncoded(i)] = biterr(tx_data(1:L), rx_bits(1:L));

    % =======================================================
    % SCHEME B: HAMMING (7,4)
    % =======================================================
    n_data_ham = floor(grid_capacity_bits / 7) * 4;
    tx_data_ham = randi([0 1], n_data_ham, 1);
    
    tx_enc = nbiot_channel_coding(tx_data_ham, 'Hamming74', 'ENCODE');
    tx_padded = [tx_enc; zeros(grid_capacity_bits - length(tx_enc), 1)];
    
    [grid, ~] = nbiot_std_mapper(tx_padded, N_sym, CellID);
    tx_wave = nbiot_ofdm_modulator(grid, N_fft, N_cp);
    
    rng(MASTER_SEED + i*100); % FREEZE CHANNEL
    rx_wave = nbiot_mcm_channel(tx_wave, N_fft, N_cp, snr);
    
    [~, rx_grid] = nbiot_rx_simple(rx_wave, N_fft, N_cp, N_sym, CellID);
    [~, H_est] = nbiot_equalizer(rx_grid, CellID);
    rx_val = rx_grid(~p_mask) ./ H_est(~p_mask);
    rx_raw = qamdemod(rx_val, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    
    rx_dec = nbiot_channel_coding(rx_raw(1:length(tx_enc)), 'Hamming74', 'DECODE');
    L = min(length(tx_data_ham), length(rx_dec));
    [~, BER_Hamming(i)] = biterr(tx_data_ham(1:L), rx_dec(1:L));

    % =======================================================
    % SCHEME C: REPETITION (3)
    % =======================================================
    n_data_rep = floor(grid_capacity_bits / 3);
    tx_data_rep = randi([0 1], n_data_rep, 1);
    
    tx_enc = nbiot_channel_coding(tx_data_rep, 'Repetition3', 'ENCODE');
    tx_padded = [tx_enc; zeros(grid_capacity_bits - length(tx_enc), 1)];
    
    [grid, ~] = nbiot_std_mapper(tx_padded, N_sym, CellID);
    tx_wave = nbiot_ofdm_modulator(grid, N_fft, N_cp);
    
    rng(MASTER_SEED + i*100); % FREEZE CHANNEL
    rx_wave = nbiot_mcm_channel(tx_wave, N_fft, N_cp, snr);
    
    [~, rx_grid] = nbiot_rx_simple(rx_wave, N_fft, N_cp, N_sym, CellID);
    [~, H_est] = nbiot_equalizer(rx_grid, CellID);
    rx_val = rx_grid(~p_mask) ./ H_est(~p_mask);
    rx_raw = qamdemod(rx_val, M, 'OutputType', 'bit', 'UnitAveragePower', true);
    
    rx_dec = nbiot_channel_coding(rx_raw(1:length(tx_enc)), 'Repetition3', 'DECODE');
    L = min(length(tx_data_rep), length(rx_dec));
    [~, BER_Repetit(i)] = biterr(tx_data_rep(1:L), rx_dec(1:L));
end

%% 3. Metrics & Visualization
schemes = {'Uncoded', 'Hamming (7,4)', 'Repetition (3)'};
rates   = [1.0,      4/7,             1/3];

% Calculate Speed
raw_kbps = (grid_capacity_bits / sim_duration_sec) / 1000;
throughputs = raw_kbps * rates;
overhead = ((1./rates) - 1) * 100;

figure('Name', 'Coding Comparison', 'Position', [100 100 900 500]);

% Plot 1: BER
subplot(1, 2, 1);
semilogy(SNR_vec, BER_Uncoded, 'k-o', 'LineWidth', 2); hold on;
semilogy(SNR_vec, BER_Hamming, 'b-s', 'LineWidth', 2);
semilogy(SNR_vec, BER_Repetit, 'r-^', 'LineWidth', 2);
grid on;
legend(schemes);
xlabel('SNR (dB)'); ylabel('BER');
title('Reliability (Lower is Better)');
ylim([1e-6 1]);

% Plot 2: Throughput
subplot(1, 2, 2);
b = bar(categorical(schemes), throughputs);
b.FaceColor = 'flat';
b.CData(1,:) = [0.2 0.2 0.2]; 
b.CData(2,:) = [0 0.4470 0.7410]; 
b.CData(3,:) = [0.6350 0.0780 0.1840]; 
ylabel('Throughput (kbps)');
title('Speed (Higher is Better)');
grid on;

%% 4. Print Table
fprintf('\n=================================================================\n');
fprintf('                 PERFORMANCE METRICS SUMMARY                     \n');
fprintf('=================================================================\n');
fprintf('%-16s | %-10s | %-12s | %-12s\n', 'Scheme', 'Code Rate', 'Overhead', 'Speed (kbps)');
fprintf('-----------------------------------------------------------------\n');
for k = 1:3
    fprintf('%-16s | %.3f      | %5.1f%%       | %8.2f\n', ...
        schemes{k}, rates(k), overhead(k), throughputs(k));
end
fprintf('-----------------------------------------------------------------\n');