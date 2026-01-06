function [rx_bits, grid_eq_1, grid_eq_2] = nbiot_dl_mimo_spatial_rx(rx_w1, rx_w2, N_fft, N_cp, N_sym, CellID, snr_db, eq_method)
% NBIOT_DL_MIMO_SPATIAL_RX Downlink Receiver with selectable Equalizer
% INPUTS:
%   ...
%   eq_method : 'MMSE', 'ZF', or 'NONE'

    %% 1. OFDM Demodulation
    [~, r1_grid] = nbiot_rx_simple(rx_w1, N_fft, N_cp, N_sym, CellID);
    [~, r2_grid] = nbiot_rx_simple(rx_w2, N_fft, N_cp, N_sym, CellID);
    
    %% 2. Channel Estimation
    mask_p0 = false(12, N_sym); mask_p1 = false(12, N_sym);
    v_shift = mod(CellID, 6);
    
    for s = 1:N_sym
        if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
            idx0 = [1, 7] + v_shift; idx0(idx0>12)=idx0(idx0>12)-12;
            mask_p0(idx0, s) = true;
            idx1 = [4, 10] + v_shift; idx1(idx1>12)=idx1(idx1>12)-12;
            mask_p1(idx1, s) = true;
        end
    end
    known_p = complex(1/sqrt(2), 1/sqrt(2));
    
    % Estimate 4 paths
    H11 = mmse_2d_block_estimate(r1_grid, mask_p0, known_p, snr_db);
    H12 = mmse_2d_block_estimate(r1_grid, mask_p1, known_p, snr_db);
    H21 = mmse_2d_block_estimate(r2_grid, mask_p0, known_p, snr_db);
    H22 = mmse_2d_block_estimate(r2_grid, mask_p1, known_p, snr_db);
    
    %% 3. MIMO Detection (Selectable)
    grid_eq_1 = complex(nan(12, N_sym), nan(12, N_sym));
    grid_eq_2 = complex(nan(12, N_sym), nan(12, N_sym));
    
    data_mask = ~(mask_p0 | mask_p1);
    noise_var = 10^(-snr_db/10); 
    
    for s = 1:N_sym
        for f = 1:12
            if data_mask(f, s)
                H = [H11(f,s), H12(f,s); H21(f,s), H22(f,s)];
                y = [r1_grid(f,s); r2_grid(f,s)];
                
                % --- EQUALIZER SELECTION ---
                switch upper(eq_method)
                    case 'MMSE'
                        % Optimal balance of noise & interference
                        W = (H' * H + noise_var * eye(2)) \ H';
                        
                    case 'ZF'
                        % Zero Forcing: Inverts channel completely
                        % W = inv(H) or pinv(H)
                        W = pinv(H);
                        
                    case 'NONE'
                        % No processing (Identity matrix)
                        % Interference will dominate
                        W = eye(2);
                        
                    otherwise
                        error('Unknown Equalizer Method');
                end
                
                x_hat = W * y;
                
                grid_eq_1(f, s) = x_hat(1);
                grid_eq_2(f, s) = x_hat(2);
            end
        end
    end
    
    %% 4. Demodulate
    rx_bits = [qamdemod(grid_eq_1(data_mask), 4, 'OutputType', 'bit', 'UnitAveragePower', true);
               qamdemod(grid_eq_2(data_mask), 4, 'OutputType', 'bit', 'UnitAveragePower', true)];
end

% (Giữ nguyên hàm con mmse_2d_block_estimate phía dưới của file cũ, copy vào đây nếu cần)
% Lưu ý: Copy phần hàm phụ mmse_2d_block_estimate từ file cũ của bạn vào cuối file này nhé.
function H_full = mmse_2d_block_estimate(grid, mask, p_val, snr_db)
    [N_subc, N_sym] = size(grid);
    H_full = zeros(N_subc, N_sym);
    block_size = 14; num_blocks = ceil(N_sym / block_size);
    safe_snr = min(snr_db, 30); beta = max(1/(10^(safe_snr/10)), 1e-4);
    tau_rms = 1e-6; fd = 70; df = 15e3; dt = 1e-3/14;
    
    for b = 1:num_blocks
        cols = (b-1)*block_size + 1 : min(b*block_size, N_sym);
        g_chunk = grid(:, cols); m_chunk = mask(:, cols);
        [pr, pc] = find(m_chunk);
        if isempty(pr), H_full(:, cols) = 1; continue; end
        
        h_ls = g_chunk(m_chunk) ./ p_val;
        dP_f = abs(pr - pr'); dP_t = abs(pc - pc');
        R_pp = (1./(1+1j*2*pi*tau_rms*dP_f*df) .* besselj(0, 2*pi*fd*dP_t*dt)) + beta*eye(length(h_ls));
        W_ls = inv(R_pp);
        
        H_blk = zeros(N_subc, length(cols));
        for si = 1:length(cols)
            dt_v = abs(si - pc')*dt;
            for fi = 1:N_subc
                 df_v = abs(fi - pr')*df;
                 r_gp = (1./(1+1j*2*pi*tau_rms*df_v)).*besselj(0,2*pi*fd*dt_v);
                 H_blk(fi, si) = r_gp * W_ls * h_ls;
            end
        end
        H_full(:, cols) = H_blk;
    end
end