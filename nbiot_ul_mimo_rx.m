function [rx_bits, rx_grid_mrc] = nbiot_ul_mimo_rx(rx_wave1, rx_wave2, N_fft, N_cp, N_sym, M, snr_db)
% NBIOT_UL_MIMO_RX Uplink Receiver using Adaptive MMSE 2D Estimation
% Added input: snr_db

    N_subc = 12;
    
    %% 1. FFT
    [grid_freq_1] = get_ul_freq_grid(rx_wave1, N_fft, N_cp, N_sym);
    [grid_freq_2] = get_ul_freq_grid(rx_wave2, N_fft, N_cp, N_sym);
    
    %% 2. Channel Estimation (MMSE 2D)
    pilot_mask = false(12, N_sym);
    for s = 1:N_sym
        if mod(s-1, 7)==3 
            pilot_mask(:, s) = true; 
        end
    end
    
    P_val = complex(1/sqrt(2), 1/sqrt(2));
    
    % Pass TRUE SNR
    H1_est = mmse_2d_block_estimate(grid_freq_1, pilot_mask, P_val, snr_db);
    H2_est = mmse_2d_block_estimate(grid_freq_2, pilot_mask, P_val, snr_db);
    
    %% 3. Maximum Ratio Combining (MRC)
    numerator = conj(H1_est) .* grid_freq_1 + conj(H2_est) .* grid_freq_2;
    denominator = abs(H1_est).^2 + abs(H2_est).^2 + 1e-10;
    
    grid_freq_combined = numerator ./ denominator;
    
    %% 4. IDFT & Demodulation
    % OLD: rx_grid_mrc = zeros(12, N_sym);
    % NEW: Initialize with NaN
    rx_grid_mrc = complex(nan(12, N_sym), nan(12, N_sym));
    
    for s = 1:N_sym
        % Check if this is a pilot symbol (Column 4 of the slot)
        % Note: pilot_mask was created in Step 2 of this file
        if any(pilot_mask(:, s))
            continue; % SKIP PILOTS! Leave them as NaN.
        end
        
        % Only perform IFFT for Data symbols
        rx_grid_mrc(:, s) = ifft(grid_freq_combined(:, s));
    end
    
    % Mask Pilots for bit extraction
    data_syms = rx_grid_mrc(~pilot_mask);
    
    % Now we can remove NaNs safely (if any slipped through) before demod
    % (Though the mask logic above should prevent them being selected)
    data_syms(isnan(data_syms)) = []; 
    
    rx_bits = qamdemod(data_syms, M, 'OutputType', 'bit', 'UnitAveragePower', true);
end

function grid = get_ul_freq_grid(wave, N_fft, N_cp, N_sym)
    sym_len = N_fft+N_cp;
    mat = reshape(wave(1:sym_len*N_sym), sym_len, N_sym);
    mat_nocp = mat(N_cp+1:end, :);
    f_grid = fftshift(fft(mat_nocp, N_fft)/sqrt(N_fft), 1);
    idx = floor((N_fft-12)/2)+1;
    grid = f_grid(idx:idx+11, :);
end

% --- HELPER: ADAPTIVE BLOCK-BASED MMSE ---
function H_full = mmse_2d_block_estimate(grid, mask, p_val, snr_db)
    [N_subc, N_sym] = size(grid);
    H_full = zeros(N_subc, N_sym);
    
    block_size = 14; 
    num_blocks = ceil(N_sym / block_size);
    
    % Update BETA based on true SNR
    beta = 1 / (10^(snr_db/10)); 
    
    tau_rms = 1e-6; fd = 91; df = 15e3; dt = 1e-3/14;
    
    for b = 1:num_blocks
        s_start = (b-1)*block_size + 1;
        s_end   = min(b*block_size, N_sym);
        cols = s_start:s_end;
        
        grid_chunk = grid(:, cols);
        mask_chunk = mask(:, cols);
        
        [p_rows, p_cols_local] = find(mask_chunk);
        
        if isempty(p_rows), H_full(:, cols) = 1; continue; end
        
        h_ls = grid_chunk(mask_chunk) ./ p_val;
        
        dP_f = abs(p_rows - p_rows');
        dP_t = abs(p_cols_local - p_cols_local');
        R_f = 1 ./ (1 + 1j * 2 * pi * tau_rms * dP_f * df); 
        R_t = besselj(0, 2 * pi * fd * dP_t * dt);          
        R_pp = (abs(R_f) .* R_t) + beta * eye(length(h_ls)); 
        W_ls = inv(R_pp); 
        
        H_block = zeros(N_subc, length(cols));
        for s_idx = 1:length(cols)
            dt_vec = abs(s_idx - p_cols_local') * dt;
            R_t_vec = besselj(0, 2 * pi * fd * dt_vec);
            for f = 1:N_subc
                df_vec = abs(f - p_rows') * df;
                R_f_vec = 1 ./ (1 + 1j * 2 * pi * tau_rms * df_vec);
                H_block(f, s_idx) = (abs(R_f_vec) .* R_t_vec) * W_ls * h_ls;
            end
        end
        H_full(:, cols) = H_block;
    end
end