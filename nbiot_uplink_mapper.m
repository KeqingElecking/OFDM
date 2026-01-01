function [resource_grid, pilot_mask] = nbiot_uplink_mapper(bits, N_sym, M)
% NBIOT_UPLINK_MAPPER Maps bits for NPUSCH Format 1 (Uplink)
% Pilots are placed at the 4th symbol of every slot.

    N_subc = 12; 
    resource_grid = zeros(N_subc, N_sym);
    pilot_mask = false(N_subc, N_sym);
    
    %% 1. Define Pilot Pattern
    % DMRS is in the 4th symbol of every slot (Indices 4, 11, 18, 25...)
    pilot_symbols = [];
    for s = 1:N_sym
        if mod(s-1, 7) == 3 
            pilot_symbols = [pilot_symbols, s];
        end
    end
    pilot_mask(:, pilot_symbols) = true;
    
    %% 2. Generate SC-FDMA Pilots (The Fix)
    % We want the pilots to be constant magnitude in the FREQUENCY domain
    % so channel estimation works on all subcarriers.
    % We define P_freq = 1, and convert to time domain for the grid.
    
    % Pilot Value (Frequency Domain)
    p_val_freq = complex(1/sqrt(2), 1/sqrt(2));
    P_vec_freq = repmat(p_val_freq, N_subc, 1);
    
    % Pilot Sequence (Time Domain)
    % We use IFFT so that when the Modulator does FFT, we get back P_vec_freq
    p_vec_time = ifft(P_vec_freq); 
    
    % Fill the pilot columns with this time-domain sequence
    for col = pilot_symbols
        resource_grid(:, col) = p_vec_time;
    end
    
    %% 3. Insert Data
    num_data_res = sum(~pilot_mask(:));
    k_bits = log2(M);
    max_bits = num_data_res * k_bits;
    
    if length(bits) > max_bits
        bits = bits(1:max_bits);
    elseif length(bits) < max_bits
        bits = [bits; randi([0 1], max_bits - length(bits), 1)];
    end
    
    data_syms = qammod(bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    resource_grid(~pilot_mask) = data_syms;
end