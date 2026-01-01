function [tx_waveform] = nbiot_scfdma_modulator(resource_grid, N_fft, N_cp)
% NBIOT_SCFDMA_MODULATOR Uplink Tx with DFT-Spreading
    
    [N_subc, N_sym] = size(resource_grid);
    
    %% 1. DFT-Spreading (The SC-FDMA Step)
    % We take the FFT of each column (12 subcarriers)
    % This moves the symbols from Time Domain representation (virtual)
    % to Frequency Domain representation.
    spread_grid = zeros(size(resource_grid));
    for s = 1:N_sym
        spread_grid(:, s) = fft(resource_grid(:, s));
    end
    
    %% 2. Map to Subcarriers (Frequency Domain)
    full_grid_shifted = zeros(N_fft, N_sym);
    
    % Map 12 subcarriers to center
    zeros_left = floor((N_fft - N_subc)/2);
    idx_start = zeros_left + 1;
    full_grid_shifted(idx_start : idx_start+N_subc-1, :) = spread_grid;
    
    %% 3. IFFT & CP (Standard OFDM generation)
    full_grid = ifftshift(full_grid_shifted, 1);
    tx_time_grid = ifft(full_grid, N_fft) * sqrt(N_fft);
    
    tx_waveform = zeros((N_fft + N_cp) * N_sym, 1);
    
    % Serial Concatenation with CP
    for s = 1:N_sym
        sym = tx_time_grid(:, s);
        cp = sym(end-N_cp+1:end);
        idx_start = (s-1)*(N_fft+N_cp) + 1;
        tx_waveform(idx_start : idx_start+N_fft+N_cp-1) = [cp; sym];
    end
end