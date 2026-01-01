function [rx_bits, rx_syms_eq] = nbiot_uplink_rx(rx_waveform, N_fft, N_cp, N_sym, M)
% NBIOT_UPLINK_RX SC-FDMA Receiver with Equalization
    
    N_subc = 12;
    
    %% 1. CP Removal & FFT
    sym_len = N_fft + N_cp;
    rx_mat = reshape(rx_waveform(1:sym_len*N_sym), sym_len, N_sym);
    rx_mat_nocp = rx_mat(N_cp+1:end, :);
    rx_fft = fft(rx_mat_nocp, N_fft) / sqrt(N_fft);
    rx_fft = fftshift(rx_fft, 1);
    
    %% 2. Extract Subcarriers
    zeros_left = floor((N_fft - N_subc)/2);
    rx_grid_freq = rx_fft(zeros_left+1 : zeros_left+N_subc, :);
    
    %% 3. Channel Estimation
    pilot_cols = [];
    for s = 1:N_sym, if mod(s-1, 7)==3, pilot_cols=[pilot_cols, s]; end; end
    
    % Reference Pilot (Frequency Domain)
    % This must match what we put in the Mapper
    p_val_freq = complex(1/sqrt(2), 1/sqrt(2));
    P_tx_freq = repmat(p_val_freq, N_subc, 1); % Vector of all non-zeros
    
    % LS Estimate at pilot columns
    H_est_cols = zeros(N_subc, length(pilot_cols));
    for i = 1:length(pilot_cols)
        col_idx = pilot_cols(i);
        Y_rx = rx_grid_freq(:, col_idx);
        
        % Division (No longer divides by zero!)
        H_est_cols(:, i) = Y_rx ./ P_tx_freq;
    end
    
    % Interpolate H across time
    H_full = zeros(N_subc, N_sym);
    
    % Handle case with only 1 pilot column (interpolation requires >=2 points)
    if length(pilot_cols) < 2
         % If only 1 pilot (e.g. N_sym=7), just replicate estimate
         H_full = repmat(H_est_cols(:,1), 1, N_sym);
    else
        for sc = 1:N_subc
            H_full(sc, :) = interp1(pilot_cols, H_est_cols(sc, :), 1:N_sym, 'linear', 'extrap');
        end
    end
    
    %% 4. Frequency Domain Equalization (FDE)
    rx_eq_freq = rx_grid_freq ./ H_full;
    
    %% 5. IDFT Despreading (Recover Time Domain Symbols)
    rx_syms_eq = zeros(N_subc, N_sym);
    for s = 1:N_sym
        rx_syms_eq(:, s) = ifft(rx_eq_freq(:, s));
    end
    
    %% 6. Demodulate
    pilot_mask = false(N_subc, N_sym);
    pilot_mask(:, pilot_cols) = true;
    
    data_syms = rx_syms_eq(~pilot_mask);
    
    % Check for NaNs just in case
    if any(isnan(data_syms))
        warning('NaNs detected in equalized symbols. Replacing with zeros.');
        data_syms(isnan(data_syms)) = 0;
    end
    
    rx_bits = qamdemod(data_syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);

end