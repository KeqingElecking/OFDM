function [rx_wave_ant1, rx_wave_ant2, H_struct] = nbiot_mimo_channel(tx_wave_ant1, tx_wave_ant2, N_fft, N_cp, snr_db)
% NBIOT_MIMO_CHANNEL Simulates a 2x2 MIMO Fading Channel
%
% Inputs:
%   tx_wave_ant1 : Signal from Tx Antenna 1 (Column Vector)
%   tx_wave_ant2 : Signal from Tx Antenna 2 (Column Vector)
%                  (For Uplink 1x2, pass zeros here)
%   snr_db       : Signal-to-Noise Ratio per receive antenna
%
% Outputs:
%   rx_wave_ant1 : Signal received at Rx Antenna 1
%   rx_wave_ant2 : Signal received at Rx Antenna 2
%   H_struct     : Struct containing channel histories (H11, H12, H21, H22)

    % 1. Setup
    % If tx_wave_ant2 is empty or zero (Uplink case), handle it
    if isempty(tx_wave_ant2) || all(tx_wave_ant2 == 0)
        tx_wave_ant2 = zeros(size(tx_wave_ant1));
        IS_UPLINK_1x2 = true;
    else
        IS_UPLINK_1x2 = false;
    end
    
    % 2. Generate 4 Independent Fading Paths
    % We use 'nbiot_mcm_channel' logic but force different random seeds internally
    % by calling it multiple times.
    % Note: We use a high SNR (100dB) in the sub-calls to get "clean" faded signals,
    % then add the actual noise at the end.
    
    % Path Tx1 -> Rx1 (h11)
    [y11, h11] = nbiot_mcm_channel(tx_wave_ant1, N_fft, N_cp, 200);
    
    % Path Tx1 -> Rx2 (h21)
    [y21, h21] = nbiot_mcm_channel(tx_wave_ant1, N_fft, N_cp, 200);
    
    if IS_UPLINK_1x2
        % Only 1 Tx antenna active
        y12 = zeros(size(y11)); h12 = zeros(size(h11));
        y22 = zeros(size(y21)); h22 = zeros(size(h21));
    else
        % Path Tx2 -> Rx1 (h12)
        [y12, h12] = nbiot_mcm_channel(tx_wave_ant2, N_fft, N_cp, 200);
        
        % Path Tx2 -> Rx2 (h22)
        [y22, h22] = nbiot_mcm_channel(tx_wave_ant2, N_fft, N_cp, 200);
    end
    
    % 3. Superposition at Receiver (Signal Combining)
    rx_clean_1 = y11 + y12; % Rx1 sees sum of Tx1 and Tx2
    rx_clean_2 = y21 + y22; % Rx2 sees sum of Tx1 and Tx2
    
    % 4. Add Noise (AWGN)
    rx_wave_ant1 = awgn(rx_clean_1, snr_db, 'measured');
    rx_wave_ant2 = awgn(rx_clean_2, snr_db, 'measured');
    
    % 5. Save Channel State (for perfect estimation debugging)
    H_struct.H11 = h11; H_struct.H12 = h12;
    H_struct.H21 = h21; H_struct.H22 = h22;

end