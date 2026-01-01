function [grid_ant0, grid_ant1] = nbiot_mimo_mapper(bits, N_sym, CellID)
% NBIOT_MIMO_MAPPER Maps bits to 2 Antenna Grids using SFBC (Transmit Diversity)

    %% 1. Setup Grids
    N_subc = 12;
    grid_ant0 = zeros(N_subc, N_sym);
    grid_ant1 = zeros(N_subc, N_sym);
    
    reserved_mask = false(N_subc, N_sym);
    
    %% 2. Insert Pilots (FIXED: Dynamic for any N_sym)
    % Standard NRS locations: Symbols 6, 7 of every 7-symbol slot.
    % Dynamic calculation allows simulation of >14 symbols.
    pilot_syms = [];
    for s = 1:N_sym
        if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
            pilot_syms = [pilot_syms, s];
        end
    end
    
    v_shift = mod(CellID, 6);
    
    for s = 1:N_sym
        if ismember(s, pilot_syms)
            % Base indices
            idx_p0 = [1, 7] + v_shift;
            idx_p1 = [4, 10] + v_shift;
            
            % Wrap
            idx_p0(idx_p0 > 12) = idx_p0(idx_p0 > 12) - 12;
            idx_p1(idx_p1 > 12) = idx_p1(idx_p1 > 12) - 12;
            
            p_val = complex(1/sqrt(2), 1/sqrt(2));
            
            % Antenna 0
            grid_ant0(idx_p0, s) = p_val;
            grid_ant0(idx_p1, s) = 0; % Muting
            
            % Antenna 1
            grid_ant1(idx_p1, s) = p_val;
            grid_ant1(idx_p0, s) = 0; % Muting
            
            % Mark Reserved
            reserved_mask(idx_p0, s) = true;
            reserved_mask(idx_p1, s) = true;
        end
    end
    
    %% 3. SFBC Data Encoding
    % FIX 2: Calculate EXACT usable pairs (accounting for skipped orphans)
    % We must mimic the loop logic to count capacity correctly.
    valid_pairs_count = 0;
    for s = 1:N_sym
        for f = 1:2:11
            % Strict check: Both f and f+1 must be free
            if ~reserved_mask(f, s) && ~reserved_mask(f+1, s)
                valid_pairs_count = valid_pairs_count + 1;
            end
        end
    end
    
    M = 4; k = 2;
    % Total bits is based on VALID pairs only
    total_bits_needed = valid_pairs_count * 2 * k; 
    
    % Handle Bits (Pad/Truncate)
    if length(bits) < total_bits_needed
        bits = [bits; randi([0 1], total_bits_needed - length(bits), 1)];
    else
        bits = bits(1:total_bits_needed);
    end
    
    qam_syms = qammod(bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    %% 4. Map to Grid
    sym_idx = 1;
    
    for s = 1:N_sym
        for f = 1:2:11
            if ~reserved_mask(f, s) && ~reserved_mask(f+1, s)
                % Safety check to prevent index error
                if sym_idx < length(qam_syms)
                    s1 = qam_syms(sym_idx);
                    s2 = qam_syms(sym_idx+1);
                    
                    % Ant0: [s1, s2]
                    grid_ant0(f, s)   = s1;
                    grid_ant0(f+1, s) = s2;
                    
                    % Ant1: [-s2*, s1*]
                    grid_ant1(f, s)   = -conj(s2);
                    grid_ant1(f+1, s) = conj(s1);
                    
                    sym_idx = sym_idx + 2;
                end
            end
        end
    end
end