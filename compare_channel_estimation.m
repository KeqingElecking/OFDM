% main_compare_equalizers.m
% COMPARISON: MMSE vs ZF vs NONE (No Filter) in NB-IoT MIMO 2x2
% ARCHITECTURE: Common Channel/Noise realization for all methods
clear; clc; close all;

%% 1. Configuration
N_fft = 64; N_cp = 32; 
N_sym = 2000; % Số lượng symbol
CellID = 2; M = 4;
SNR_vec = 0:4:32; % Bước nhảy 4dB

Methods = {'MMSE', 'ZF', 'NONE'};
Num_Methods = length(Methods);
Num_SNR = length(SNR_vec);

% Lưu kết quả: [Hàng = Method, Cột = SNR]
BER_DL_Results = zeros(Num_Methods, Num_SNR);
BER_UL_Results = zeros(Num_Methods, Num_SNR);

MASTER_SEED = 12345; 

fprintf('================================================================\n');
fprintf('  MIMO EQUALIZER COMPARISON (MMSE vs ZF vs NONE)  \n');
fprintf('================================================================\n');
fprintf('%-10s | %-10s | %-12s | %-12s\n', 'SNR (dB)', 'Method', 'BER DL', 'BER UL');
fprintf('----------------------------------------------------------------\n');

%% 2. Simulation Loop (Loop SNR first -> Generate Data -> Loop Methods)
for s_idx = 1:Num_SNR
    snr = SNR_vec(s_idx);
    
    % ---------------------------------------------------------
    % STEP A: FREEZE CHANNEL & GENERATE COMMON DATA/WAVEFORMS
    % ---------------------------------------------------------
    % Reset seed tại mỗi mức SNR để đảm bảo tái lập được kết quả
    rng(MASTER_SEED + s_idx); 
    
    % --- 1. COMMON DOWNLINK GENERATION ---
    total_dl_bits = floor((12*N_sym - 2*N_sym/7*2) * 2 * 2); 
    tx_bits_dl = randi([0 1], total_dl_bits, 1);
    
    [g_dl_0, g_dl_1] = nbiot_dl_spatial_mapper(tx_bits_dl, N_sym, CellID);
    tx_dl_0 = nbiot_ofdm_modulator(g_dl_0, N_fft, N_cp);
    tx_dl_1 = nbiot_ofdm_modulator(g_dl_1, N_fft, N_cp);
    
    % Truyền qua kênh (Tạo ra rx_dl chung cho mọi thuật toán)
    [rx_dl_0, rx_dl_1, ~] = nbiot_mimo_channel(tx_dl_0, tx_dl_1, N_fft, N_cp, snr);
    
    % --- 2. COMMON UPLINK GENERATION ---
    total_ul_bits = floor((12*N_sym - 12*N_sym/7) * 2 * 2); 
    tx_bits_ul = randi([0 1], total_ul_bits, 1);
    
    [g_ul_0, g_ul_1] = nbiot_ul_spatial_mapper(tx_bits_ul, N_sym, M);
    tx_ul_0 = nbiot_scfdma_modulator(g_ul_0, N_fft, N_cp);
    tx_ul_1 = nbiot_scfdma_modulator(g_ul_1, N_fft, N_cp);
    
    % Truyền qua kênh (Tạo ra rx_ul chung cho mọi thuật toán)
    [rx_ul_0, rx_ul_1, ~] = nbiot_mimo_channel(tx_ul_0, tx_ul_1, N_fft, N_cp, snr);
    
    % ---------------------------------------------------------
    % STEP B: TEST EACH METHOD ON THE SAME RECEIVED SIGNAL
    % ---------------------------------------------------------
    for m_idx = 1:Num_Methods
        current_method = Methods{m_idx};
        
        % --- Downlink Processing ---
        [rx_bits_dl, ~, ~] = nbiot_dl_mimo_spatial_rx(...
            rx_dl_0, rx_dl_1, N_fft, N_cp, N_sym, CellID, snr, current_method);
        
        L_dl = min(length(tx_bits_dl), length(rx_bits_dl));
        [~, ber_dl] = biterr(tx_bits_dl(1:L_dl), rx_bits_dl(1:L_dl));
        BER_DL_Results(m_idx, s_idx) = ber_dl;
        
        % --- Uplink Processing ---
        [rx_bits_ul, ~, ~] = nbiot_ul_mimo_spatial_rx(...
            rx_ul_0, rx_ul_1, N_fft, N_cp, N_sym, M, snr, current_method);
        
        L_ul = min(length(tx_bits_ul), length(rx_bits_ul));
        [~, ber_ul] = biterr(tx_bits_ul(1:L_ul), rx_bits_ul(1:L_ul));
        BER_UL_Results(m_idx, s_idx) = ber_ul;
        
        % Log kết quả ngay lập tức
        fprintf('%-10d | %-10s | %-12.2e | %-12.2e\n', snr, current_method, ber_dl, ber_ul);
    end
    fprintf('----------------------------------------------------------------\n');
end

%% 3. Visualization
figure('Name', 'MIMO Equalizer Comparison', 'Position', [100 100 1200 500]);

% Plot 1: Downlink Results
subplot(1, 2, 1);
semilogy(SNR_vec, BER_DL_Results(1, :), 'b-o', 'LineWidth', 2, 'DisplayName', 'MMSE (Optimal)'); hold on;
semilogy(SNR_vec, BER_DL_Results(2, :), 'r-s', 'LineWidth', 2, 'DisplayName', 'Zero Forcing (ZF)');
semilogy(SNR_vec, BER_DL_Results(3, :), 'k--^', 'LineWidth', 2, 'DisplayName', 'None (No Filter)');
grid on; xlabel('SNR (dB)'); ylabel('BER'); title('Downlink (OFDMA) Performance');
legend('Location', 'SouthWest'); ylim([1e-5 1]);

% Plot 2: Uplink Results
subplot(1, 2, 2);
semilogy(SNR_vec, BER_UL_Results(1, :), 'b-o', 'LineWidth', 2, 'DisplayName', 'MMSE (Optimal)'); hold on;
semilogy(SNR_vec, BER_UL_Results(2, :), 'r-s', 'LineWidth', 2, 'DisplayName', 'Zero Forcing (ZF)');
semilogy(SNR_vec, BER_UL_Results(3, :), 'k--^', 'LineWidth', 2, 'DisplayName', 'None (No Filter)');
grid on; xlabel('SNR (dB)'); ylabel('BER'); title('Uplink (SC-FDMA) Performance');
legend('Location', 'SouthWest'); ylim([1e-5 1]);

fprintf('\nSimulation Complete.\n');