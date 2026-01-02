function [rx_bits, rx_grid] = nbiot_rx_simple(rx_waveform, N_fft, N_cp, N_sym, CellID)
% NBIOT_RX_SIMPLE Basic Receiver (CP Removal -> FFT -> Demodulation)
%
% FIX APPLIED: Dynamic Pilot Mask generation for any duration (N_sym)

    %% 1. Serial to Parallel & CP Removal
    sym_len = N_fft + N_cp;
    
    % Ensure length matches expected symbols
    if length(rx_waveform) < sym_len * N_sym
        error('Input waveform is too short for N_sym symbols.');
    end
    
    % Reshape to [Symbol_Length x Num_Symbols]
    rx_matrix_cp = reshape(rx_waveform(1:sym_len*N_sym), sym_len, N_sym);
    
    % Remove CP (Discard first N_cp samples of each column)
    rx_matrix_time = rx_matrix_cp(N_cp+1:end, :);
    
    %% 2. FFT & Shift
    % FFT down columns, normalize by sqrt(N_fft)
    rx_matrix_fft = fft(rx_matrix_time, N_fft) / sqrt(N_fft);
    rx_grid_full = fftshift(rx_matrix_fft, 1);
    
    %% 3. Extract NB-IoT Subcarriers
    % Modulator placed 12 subcarriers in the center
    total_zeros = N_fft - 12;
    zeros_left = floor(total_zeros / 2);
    
    idx_start = zeros_left + 1;
    idx_end   = idx_start + 11;
    
    rx_grid = rx_grid_full(idx_start:idx_end, :);
    
    %% 4. Pilot Mask Generation (FIXED: Dynamic)
    % Instead of hardcoding [6, 7, 13, 14], we check every symbol using modulo math.
    % Standard NB-IoT: Pilots are in the last 2 symbols of every 7-symbol slot.
    % (Indices 5 and 6 if counting from 0).
    
    pilot_mask = false(12, N_sym);
    v_shift = mod(CellID, 6);
    
    for s = 1:N_sym
        % Check if current symbol index 's' is a pilot symbol
        if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
            % It is a pilot symbol! Calculate positions.
            p_subs = [1, 7] + v_shift; % Base locations
            
            % Wrap around frequency domain if needed
            p_subs(p_subs > 12) = p_subs(p_subs > 12) - 12; 
            
            pilot_mask(p_subs, s) = true;
        end
    end
    
    %% 5. Demodulation (QPSK)
    % Extract ONLY Data REs (Skip Pilots)
    rx_data_syms = rx_grid(~pilot_mask);
    
    % QPSK Demodulation
    rx_bits = qamdemod(rx_data_syms, 4, 'OutputType', 'bit', 'UnitAveragePower', true);

end