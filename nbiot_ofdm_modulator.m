function [tx_waveform] = nbiot_ofdm_modulator(resource_grid, N_fft, N_cp)
% NBIOT_OFDM_MODULATOR Performs Zero Insertion, IFFT, and CP Addition
%
% Inputs:
%   resource_grid : Complex grid [12 x N_sym] (Output from nbiot_dl_modulator)
%   N_fft         : FFT Size (Standard NB-IoT simulation uses 64 or 128)
%   N_cp          : Cyclic Prefix length (in samples)
%
% Output:
%   tx_waveform   : Serial Time-Domain Vector (Column)

    %% 1. Input Validation
    [num_subc, num_sym] = size(resource_grid);
    
    if num_subc ~= 12
        error('NB-IoT Resource Grid must have exactly 12 subcarriers.');
    end
    
    if N_fft < 12
        error('FFT size must be larger than 12.');
    end

    %% 2. Zero Insertion (Spectral Mapping)
    % We map the 12 active subcarriers to the CENTER of the spectrum.
    % Formula: [Zeros_Left | Data | Zeros_Right]
    
    % Calculate number of zeros needed
    total_zeros = N_fft - 12;
    zeros_left = floor(total_zeros / 2);
    zeros_right = ceil(total_zeros / 2);
    
    % Create the Full Frequency Grid (64 x N_sym)
    % We initialize with zeros
    full_grid_shifted = zeros(N_fft, num_sym);
    
    % Insert Data in the middle
    % Indices: (Zeros_Left + 1) to (Zeros_Left + 12)
    idx_start = zeros_left + 1;
    idx_end = idx_start + 11;
    
    full_grid_shifted(idx_start:idx_end, :) = resource_grid;
    
    %% 3. IFFT Processing
    % CRITICAL: ifftshift
    % We constructed the grid visually [Left...Center...Right].
    % MATLAB's IFFT expects [DC...Pos...Neg].
    % We must swap the halves before IFFT.
    full_grid = ifftshift(full_grid_shifted, 1);
    
    % Perform IFFT down the columns
    % We scale by sqrt(N_fft) to preserve energy (optional, but standard in LTE)
    tx_time_grid = ifft(full_grid, N_fft) * sqrt(N_fft);
    
    %% 4. Guard Insertion (Cyclic Prefix)
    % We append the end of the symbol to the beginning.
    
    % Create buffer for time symbols with CP
    tx_time_with_cp = zeros(N_fft + N_cp, num_sym);
    
    for s = 1:num_sym
        % Extract current symbol
        sym_time = tx_time_grid(:, s);
        
        % Copy tail (CP)
        cp_segment = sym_time(end - N_cp + 1 : end);
        
        % Concatenate
        tx_time_with_cp(:, s) = [cp_segment; sym_time];
    end
    
    %% 5. Parallel to Serial Conversion
    % Flatten the matrix into a single column vector stream
    tx_waveform = tx_time_with_cp(:);

end