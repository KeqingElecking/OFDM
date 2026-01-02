function [resource_grid, pilot_mask] = nbiot_std_mapper(bits, N_sym, CellID)
% NBIOT_STD_MAPPER Maps bits using Standard NB-IoT Pilot Patterns (NRS)
% Updated to support N_sym > 14 (Dynamic Pilot Calculation)

    %% 1. Setup
    N_subc = 12; 
    resource_grid = zeros(N_subc, N_sym);
    pilot_mask = false(N_subc, N_sym);
    
    %% 2. Define Standard NRS Pilot Pattern (DYNAMIC)
    % Rules: Pilots appear in the last 2 symbols of each 7-symbol slot.
    % Indices (0-based in slot): 5, 6
    
    pilot_symbols = [];
    for s = 1:N_sym
        if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
            pilot_symbols = [pilot_symbols, s];
        end
    end
    
    v_shift = mod(CellID, 6);
    
    for s = 1:N_sym
        if ismember(s, pilot_symbols)
            % Determine pilot subcarriers for this symbol
            % Basic pattern: indices [1, 7] shifted by v_shift
            p_subs = [1, 7] + v_shift;
            
            % Wrap around
            p_subs(p_subs > 12) = p_subs(p_subs > 12) - 12;
            
            % Mark grid
            pilot_mask(p_subs, s) = true;
        end
    end
    
    %% 3. Insert Pilots
    p_val = complex(1/sqrt(2), 1/sqrt(2));
    resource_grid(pilot_mask) = p_val;
    
    %% 4. Insert Data
    % Calculate available Data Resource Elements (REs)
    num_pilot_res = sum(pilot_mask(:));
    num_data_res = (N_subc * N_sym) - num_pilot_res;
    
    % NB-IoT Downlink is QPSK (2 bits/symbol)
    M = 4; k = 2;
    max_bits = num_data_res * k;
    
    % Handle bit count (Pad or Truncate)
    if length(bits) > max_bits
        bits = bits(1:max_bits);
    elseif length(bits) < max_bits
        % If input is shorter, pad with random bits
        % (This happens if capacity calc was wrong in script)
        padding = randi([0 1], max_bits - length(bits), 1);
        bits = [bits; padding];
    end
    
    % Modulate Data
    data_syms = qammod(bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    % Map to Grid (Filling only non-pilot spots)
    resource_grid(~pilot_mask) = data_syms;

end