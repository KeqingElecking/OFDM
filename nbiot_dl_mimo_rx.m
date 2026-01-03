function [rx_bits, grid_eq] = nbiot_dl_mimo_rx(rx_wave, N_fft, N_cp, N_sym, CellID, snr_db)
% NBIOT_DL_MIMO_RX Downlink Receiver 
% Features: Adaptive MMSE 2D Estimation & NaN-based Pilot Removal

    %% 1. OFDM Demodulation
    [~, rx_grid] = nbiot_rx_simple(rx_wave, N_fft, N_cp, N_sym, CellID);
    
    %% 2. MIMO Channel Estimation (MMSE 2D)
    mask_p1 = false(12, N_sym);
    mask_p2 = false(12, N_sym);
    v_shift = mod(CellID, 6);
    
    for s = 1:N_sym
        if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
            idx1 = [1, 7] + v_shift;
            idx1(idx1>12) = idx1(idx1>12)-12;
            mask_p1(idx1, s) = true;
            
            idx2 = [4, 10] + v_shift;
            idx2(idx2>12) = idx2(idx2>12)-12;
            mask_p2(idx2, s) = true;
        end
    end
    
    known_p = complex(1/sqrt(2), 1/sqrt(2));
    
    % Pass TRUE SNR to the estimator
    H1_est = mmse_2d_block_estimate(rx_grid, mask_p1, known_p, snr_db);
    H2_est = mmse_2d_block_estimate(rx_grid, mask_p2, known_p, snr_db);
    
    %% 3. SFBC Decoding
    % Initialize with NaN so plots automatically ignore invalid points
    grid_eq = complex(nan(12, N_sym), nan(12, N_sym));
    
    all_pilots = mask_p1 | mask_p2;
    data_mask = false(12, N_sym); 
    
    for s = 1:N_sym
        for f = 1:2:11
            % Check pair for pilots: Process ONLY if both are data
            if ~all_pilots(f, s) && ~all_pilots(f+1, s)
                
                r1 = rx_grid(f, s);
                r2 = rx_grid(f+1, s);
                
                h1 = H1_est(f, s);
                h2 = H2_est(f, s);
                
                % --- RESTORED MATH START ---
                % Alamouti Combiner (SFBC)
                s1_hat = conj(h1)*r1 + h2*conj(r2);
                s2_hat = conj(h2)*r1 - h1*conj(r2);
                % --- RESTORED MATH END ---
                
                pow_norm = abs(h1)^2 + abs(h2)^2 + 1e-10;
                
                grid_eq(f, s)   = s1_hat / pow_norm;
                grid_eq(f+1, s) = s2_hat / pow_norm;
                
                data_mask(f, s) = true;
                data_mask(f+1, s) = true;
            end
        end
    end
    
    %% 4. Demodulate
    % rx_data will be a column vector of valid QPSK symbols
    rx_data = grid_eq(data_mask);
    rx_bits = qamdemod(rx_data, 4, 'OutputType', 'bit', 'UnitAveragePower', true);
end

% --- HELPER: ADAPTIVE BLOCK-BASED MMSE ---
function H_full = mmse_2d_block_estimate(grid, mask, p_val, snr_db)
    [N_subc, N_sym] = size(grid);
    H_full = zeros(N_subc, N_sym);
    
    block_size = 14; 
    num_blocks = ceil(N_sym / block_size);
    
    % Stability Fix: Cap SNR to prevent matrix singularity
    safe_snr = min(snr_db, 30);
    beta = 1 / (10^(safe_snr/10)); 
    beta = max(beta, 1e-4); % Minimum floor
    
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