% run_full_spatial_mimo.m
% Comprehensive NB-IoT TRUE MIMO (Spatial Multiplexing) Simulation
% Features: 2x2 Spatial Mux, MMSE Receiver, Channel Physics, REAL-TIME LOGGING
clear; clc; close all;

%% 1. Configuration
N_fft = 64; N_cp = 32; 
N_sym = 10000; 
CellID = 2; M = 4;
SNR_vec = 0:2:32; 

BER_DL = zeros(size(SNR_vec));
BER_UL = zeros(size(SNR_vec));
MASTER_SEED = 55555; 

disp('Running 2x2 Spatial Multiplexing Loop (Standard MMSE)...');
disp('-------------------------------------------------------------');
fprintf('%-10s | %-18s | %-18s\n', 'SNR (dB)', 'BER Downlink', 'BER Uplink');
disp('-------------------------------------------------------------');

%% 2. Simulation Loop
for i = 1:length(SNR_vec)
    snr = SNR_vec(i);
    
    % FREEZE CHANNEL PHYSICS (Critical for fair comparison)
    rng(MASTER_SEED); 
    
    % =========================================================
    % A. DOWNLINK (Spatial Mux 2x2 OFDM)
    % =========================================================
    total_dl_bits = floor((12*N_sym - 2*N_sym/7*2) * 2 * 2); 
    tx_bits_dl = randi([0 1], total_dl_bits, 1);
    
    [g_dl_0, g_dl_1] = nbiot_dl_spatial_mapper(tx_bits_dl, N_sym, CellID);
    
    tx_dl_0 = nbiot_ofdm_modulator(g_dl_0, N_fft, N_cp);
    tx_dl_1 = nbiot_ofdm_modulator(g_dl_1, N_fft, N_cp);
    
    [rx_dl_0, rx_dl_1, H_struct_dl] = nbiot_mimo_channel(tx_dl_0, tx_dl_1, N_fft, N_cp, snr);
    
    % --- UPDATED CALL: Added 'MMSE' argument ---
    [rx_bits_dl, eq_dl_1, eq_dl_2] = nbiot_dl_mimo_spatial_rx(...
        rx_dl_0, rx_dl_1, N_fft, N_cp, N_sym, CellID, snr, 'MMSE');
    
    L = min(length(tx_bits_dl), length(rx_bits_dl));
    [~, BER_DL(i)] = biterr(tx_bits_dl(1:L), rx_bits_dl(1:L));
    
    % =========================================================
    % B. UPLINK (Spatial Mux 2x2 SC-FDMA)
    % =========================================================
    total_ul_bits = floor((12*N_sym - 12*N_sym/7) * 2 * 2); 
    tx_bits_ul = randi([0 1], total_ul_bits, 1);
    
    [g_ul_0, g_ul_1] = nbiot_ul_spatial_mapper(tx_bits_ul, N_sym, M);
    
    tx_ul_0 = nbiot_scfdma_modulator(g_ul_0, N_fft, N_cp);
    tx_ul_1 = nbiot_scfdma_modulator(g_ul_1, N_fft, N_cp);
    
    [rx_ul_0, rx_ul_1, H_struct_ul] = nbiot_mimo_channel(tx_ul_0, tx_ul_1, N_fft, N_cp, snr);
    
    % --- UPDATED CALL: Added 'MMSE' argument ---
    [rx_bits_ul, eq_ul_1, eq_ul_2] = nbiot_ul_mimo_spatial_rx(...
        rx_ul_0, rx_ul_1, N_fft, N_cp, N_sym, M, snr, 'MMSE');
    
    L_ul = min(length(tx_bits_ul), length(rx_bits_ul));
    [~, BER_UL(i)] = biterr(tx_bits_ul(1:L_ul), rx_bits_ul(1:L_ul));
    
    % =========================================================
    % LOGGING RESULTS TO CONSOLE
    % =========================================================
    fprintf('%-10d | %-18.2e | %-18.2e\n', snr, BER_DL(i), BER_UL(i));

    % --- SAVE DATA FOR PLOTS (at 20dB) ---
    if snr == 20
        saved_dl_L1 = eq_dl_1; saved_dl_L2 = eq_dl_2;
        saved_ul_L1 = eq_ul_1; saved_ul_L2 = eq_ul_2;
        saved_H_dl = H_struct_dl.H11;
        saved_H_ul = H_struct_ul.H11;
        saved_wave_tx_dl = tx_dl_0; saved_wave_rx_dl = rx_dl_0;
        saved_wave_tx_ul = tx_ul_0; saved_wave_rx_ul = rx_ul_0;
    end
end
disp('-------------------------------------------------------------');
disp('Simulation Complete.');

%% 3. Visualization

% Figure 1: BER
figure('Name', 'MIMO Spatial Multiplexing Performance', 'Position', [100 100 600 400]);
semilogy(SNR_vec, BER_DL, 'b-o', 'LineWidth', 2); hold on;
semilogy(SNR_vec, BER_UL, 'r-s', 'LineWidth', 2);
grid on; xlabel('SNR (dB)'); ylabel('BER');
title('2x2 Spatial Multiplexing (True MIMO)');
legend('Downlink (OFDM)', 'Uplink (SC-FDMA)'); ylim([1e-5 1]);

% Figure 2: Constellations (20dB)
figure('Name', 'Recovered Layers (20dB)', 'Position', [750 100 800 400]);
subplot(1,2,1);
plot(saved_dl_L1(:), 'b.'); hold on;
plot([1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2), 'ro', 'LineWidth', 2);
title('DL Layer 1 (OFDM)'); axis square; grid on;

subplot(1,2,2);
plot(saved_ul_L1(:), 'm.'); hold on;
plot([1+1j, 1-1j, -1+1j, -1-1j]/sqrt(2), 'ro', 'LineWidth', 2);
title('UL Layer 1 (SC-FDMA)'); axis square; grid on;

% Figure 3: Channel Physics (Dual Plot)
figure('Name', 'Channel Fading Profile', 'Position', [400 550 600 700]);
[~, n_syms_plot] = size(saved_H_dl);
subcarrier_spacing = 15e3; 
freq_axis_khz = ((-6:5) * subcarrier_spacing) / 1000; 
center_idx = floor((N_fft - 12)/2) + 1; 
active_idx = center_idx : (center_idx + 11);

% DL Channel
H_freq_dl = zeros(N_fft, n_syms_plot);
for s = 1:n_syms_plot, H_freq_dl(:, s) = fft(saved_H_dl(:, s), N_fft); end
H_active_dl = fftshift(H_freq_dl, 1); H_active_dl = H_active_dl(active_idx, :);

subplot(2,1,1);
surf(1:n_syms_plot, freq_axis_khz, abs(H_active_dl)); 
shading interp; colorbar; view(-45, 30);
title('Downlink Channel Magnitude'); ylabel('Freq (kHz)'); zlabel('|H|'); axis tight;

% UL Channel
H_freq_ul = zeros(N_fft, n_syms_plot);
for s = 1:n_syms_plot, H_freq_ul(:, s) = fft(saved_H_ul(:, s), N_fft); end
H_active_ul = fftshift(H_freq_ul, 1); H_active_ul = H_active_ul(active_idx, :);

subplot(2,1,2);
surf(1:n_syms_plot, freq_axis_khz, abs(H_active_ul)); 
shading interp; colorbar; view(-45, 30);
title('Uplink Channel Magnitude'); xlabel('Symbol'); ylabel('Freq (kHz)'); zlabel('|H|'); axis tight;

% --- FIGURE 4: SPECTRUM ANALYZER ---
figure('Name', 'Spectrum Analyzer', 'Position', [1050 100 800 600]);
Fs = N_fft * 15e3; 
f_axis_plot = linspace(-Fs/2, Fs/2, 4096)/1000;

% Helper for PSD
calc_psd = @(x) 10*log10(abs(fftshift(fft(x(1:4096)))).^2);

subplot(2,2,1); plot(f_axis_plot, calc_psd(saved_wave_tx_dl), 'b'); grid on; title('DL Tx'); xlim([-400 400]);
subplot(2,2,2); plot(f_axis_plot, calc_psd(saved_wave_rx_dl), 'b'); grid on; title('DL Rx'); xlim([-400 400]);
subplot(2,2,3); plot(f_axis_plot, calc_psd(saved_wave_tx_ul), 'm'); grid on; title('UL Tx'); xlim([-400 400]);
subplot(2,2,4); plot(f_axis_plot, calc_psd(saved_wave_rx_ul), 'm'); grid on; title('UL Rx'); xlim([-400 400]);

%% 4. Effective Spectral Efficiency
delta_f = 15e3; N_subc = 12; B = N_subc * delta_f;
Fs = N_fft * delta_f; T_sym = (N_fft + N_cp) / Fs; Total_Time = N_sym * T_sym;

Rb_DL_bps = length(tx_bits_dl) / Total_Time; Reff_DL = Rb_DL_bps / B;
Rb_UL_bps = length(tx_bits_ul) / Total_Time; Reff_UL = Rb_UL_bps / B;

fprintf('\n=== REPORT ===\n');
fprintf('DL Efficiency: %.2f bits/s/Hz\n', Reff_DL);
fprintf('UL Efficiency: %.2f bits/s/Hz\n', Reff_UL);