function [rx_bits, rx_grid] = nbiot_rx_simple(rx_waveform, N_fft, N_cp, N_sym, CellID)
% NBIOT_RX_SIMPLE Basic Receiver (CP Removal -> FFT -> Demodulation)
% note: Does NOT perform channel equalization.
%
% Inputs:
%   rx_waveform : Serial time-domain signal
%   N_fft       : FFT Size (e.g., 64)
%   N_cp        : Cyclic Prefix Length
%   N_sym       : Number of symbols expected
%   CellID      : To identify and skip Pilot locations
%
% Outputs:
%   rx_bits     : Demodulated bits (column vector)
%   rx_grid     : 12 x N_sym complex grid (for debugging/equalizer)

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
    % FFT down columns
    % Scale by 1/sqrt(N_fft) to normalize the modulator's gain
    rx_matrix_fft = fft(rx_matrix_time, N_fft) / sqrt(N_fft);
    
    % FFT Shift (Move DC from index 1 back to center)
    rx_grid_full = fftshift(rx_matrix_fft, 1);
    
    %% 3. Extract NB-IoT Subcarriers
    % Modulator placed 12 subcarriers in the center
    total_zeros = N_fft - 12;
    zeros_left = floor(total_zeros / 2);
    
    idx_start = zeros_left + 1;
    idx_end   = idx_start + 11;
    
    % Extract the 12xN_sym grid
    rx_grid = rx_grid_full(idx_start:idx_end, :);
    
    %% 4. Pilot Mask Generation (to skip pilots)
    % Re-create the pilot mask using standard NB-IoT logic
    % (Logic copied from nbiot_std_mapper)
    pilot_mask = false(12, N_sym);
    pilot_symbols = [6, 7, 13, 14]; % Standard NRS locations
    v_shift = mod(CellID, 6);
    
    for s = 1:N_sym
        if ismember(s, pilot_symbols)
            p_subs = [1, 7] + v_shift; % Base locations
            p_subs(p_subs > 12) = p_subs(p_subs > 12) - 12; % Wrap
            pilot_mask(p_subs, s) = true;
        end
    end
    
    %% 5. Demodulation (QPSK)
    % Extract ONLY Data REs (Skip Pilots)
    rx_data_syms = rx_grid(~pilot_mask);
    
    % QPSK Demodulation (Hard Decision)
    % UnitAveragePower=true matches modulator scaling
    rx_bits = qamdemod(rx_data_syms, 4, 'OutputType', 'bit', 'UnitAveragePower', true);

end