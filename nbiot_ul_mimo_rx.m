function [rx_bits, rx_grid_mrc] = nbiot_ul_mimo_rx(rx_wave1, rx_wave2, N_fft, N_cp, N_sym, M)
% NBIOT_UL_MIMO_RX Uplink Receiver for 1x2 SIMO (MRC)
% 
    N_subc = 12;
    
    %% 1. FFT & Extract Subcarriers (Per Antenna)
    [grid_freq_1] = get_ul_freq_grid(rx_wave1, N_fft, N_cp, N_sym);
    [grid_freq_2] = get_ul_freq_grid(rx_wave2, N_fft, N_cp, N_sym);
    
    %% 2. Channel Estimation (Per Antenna)
    % Pilot Logic (Symbol 4)
    pilot_cols = [];
    for s = 1:N_sym, if mod(s-1, 7)==3, pilot_cols=[pilot_cols, s]; end; end
    
    % Uplink Pilot Value (Freq Domain)
    P_tx = repmat(complex(1/sqrt(2), 1/sqrt(2)), 12, 1);
    
    H1_est = estimate_ul_channel(grid_freq_1, pilot_cols, P_tx, N_sym);
    H2_est = estimate_ul_channel(grid_freq_2, pilot_cols, P_tx, N_sym);
    
    %% 3. Maximum Ratio Combining (MRC)
    % H_mrc = H1* Y1 + H2* Y2  / (|H1|^2 + |H2|^2)
    
    numerator = conj(H1_est) .* grid_freq_1 + conj(H2_est) .* grid_freq_2;
    denominator = abs(H1_est).^2 + abs(H2_est).^2 + 1e-10;
    
    grid_freq_combined = numerator ./ denominator;
    
    %% 4. IDFT & Demodulation
    rx_grid_mrc = zeros(12, N_sym);
    for s = 1:N_sym
        rx_grid_mrc(:, s) = ifft(grid_freq_combined(:, s));
    end
    
    % Mask Pilots
    p_mask = false(12, N_sym);
    p_mask(:, pilot_cols) = true;
    
    data_syms = rx_grid_mrc(~p_mask);
    
    % Handle potential NaNs from empty pilot slots
    data_syms(isnan(data_syms)) = 0;
    
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

function H_full = estimate_ul_channel(grid, p_cols, P_val, N_sym)
    H_full = zeros(12, N_sym);
    H_est_cols = grid(:, p_cols) ./ P_val;
    if length(p_cols) < 2
         H_full = repmat(H_est_cols(:,1), 1, N_sym);
    else
        for sc=1:12
            H_full(sc,:) = interp1(p_cols, H_est_cols(sc,:), 1:N_sym, 'linear', 'extrap');
        end
    end
end