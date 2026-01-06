function [rx_bits, grid_eq_1, grid_eq_2] = nbiot_ul_mimo_spatial_rx(rx_w1, rx_w2, N_fft, N_cp, N_sym, M, snr_db, eq_method)
% NBIOT_UL_MIMO_SPATIAL_RX Uplink Receiver with selectable Equalizer

    %% 1. FFT
    [grid_f1] = get_ul_freq_grid(rx_w1, N_fft, N_cp, N_sym);
    [grid_f2] = get_ul_freq_grid(rx_w2, N_fft, N_cp, N_sym);
    
    %% 2. Channel Estimation
    mask_p0 = false(12, N_sym); mask_p1 = false(12, N_sym);
    for s = 1:N_sym
        if mod(s-1, 7) == 3
            mask_p0(1:2:12, s) = true; 
            mask_p1(2:2:12, s) = true; 
        end
    end
    known_p = complex(1/sqrt(2), 1/sqrt(2)); 
    
    % Estimate 4 paths
    H11 = mmse_2d_block_estimate(grid_f1, mask_p0, known_p, snr_db);
    H12 = mmse_2d_block_estimate(grid_f1, mask_p1, known_p, snr_db);
    H21 = mmse_2d_block_estimate(grid_f2, mask_p0, known_p, snr_db);
    H22 = mmse_2d_block_estimate(grid_f2, mask_p1, known_p, snr_db);
    
    %% 3. MIMO Detection (Frequency Domain)
    grid_eq_f1 = zeros(12, N_sym);
    grid_eq_f2 = zeros(12, N_sym);
    
    noise_var = 10^(-snr_db/10); 
    
    for s = 1:N_sym
        for f = 1:12
             H = [H11(f,s), H12(f,s); H21(f,s), H22(f,s)];
             y = [grid_f1(f,s); grid_f2(f,s)];
             
             % --- EQUALIZER SELECTION ---
             switch upper(eq_method)
                case 'MMSE'
                    W = (H' * H + noise_var * eye(2)) \ H';
                case 'ZF'
                    W = pinv(H);
                case 'NONE'
                    W = eye(2);
                otherwise
                    error('Unknown Method');
             end
             
             x_hat = W * y;
             
             grid_eq_f1(f, s) = x_hat(1);
             grid_eq_f2(f, s) = x_hat(2);
        end
    end
    
    %% 4. IDFT (SC-FDMA Despreading)
    grid_eq_1 = complex(nan(12, N_sym), nan(12, N_sym));
    grid_eq_2 = complex(nan(12, N_sym), nan(12, N_sym));
    
    for s = 1:N_sym
        if any(mask_p0(:,s)) || any(mask_p1(:,s)), continue; end
        grid_eq_1(:, s) = ifft(grid_eq_f1(:, s));
        grid_eq_2(:, s) = ifft(grid_eq_f2(:, s));
    end
    
    %% 5. Demodulate
    syms_1 = grid_eq_1(~isnan(grid_eq_1));
    syms_2 = grid_eq_2(~isnan(grid_eq_2));
    
    rx_bits = [qamdemod(syms_1, M, 'OutputType', 'bit', 'UnitAveragePower', true);
               qamdemod(syms_2, M, 'OutputType', 'bit', 'UnitAveragePower', true)];
end

% (Copy hàm phụ get_ul_freq_grid và mmse_2d_block_estimate từ code cũ vào đây)
function grid = get_ul_freq_grid(wave, N_fft, N_cp, N_sym)
    sym_len = N_fft+N_cp;
    mat = reshape(wave(1:sym_len*N_sym), sym_len, N_sym);
    f_grid = fftshift(fft(mat(N_cp+1:end, :), N_fft)/sqrt(N_fft), 1);
    idx = floor((N_fft-12)/2)+1;
    grid = f_grid(idx:idx+11, :);
end

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