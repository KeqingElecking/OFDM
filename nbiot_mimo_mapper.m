function [grid_ant0, grid_ant1] = nbiot_mimo_mapper(bits, N_sym, CellID)
% NBIOT_MIMO_MAPPER Maps bits to 2 Antenna Grids using SFBC (Transmit Diversity)
%
% Inputs:
%   bits    : Data bits (Column vector)
%   N_sym   : Number of symbols (14)
%   CellID  : Physical Cell ID
%
% Outputs:
%   grid_ant0 : 12 x N_sym complex grid for Antenna 0
%   grid_ant1 : 12 x N_sym complex grid for Antenna 1

    %% 1. Setup Grids
    N_subc = 12;
    grid_ant0 = zeros(N_subc, N_sym);
    grid_ant1 = zeros(N_subc, N_sym);
    
    % Track which REs are occupied by pilots (to skip them for data)
    reserved_mask = false(N_subc, N_sym);
    
    %% 2. Insert Pilots (Orthogonal)
    % For simulation, we use a "Pilot on A / Null on B" approach.
    % Ant 0 Pilots: Standard NRS locations
    % Ant 1 Pilots: Standard NRS locations shifted by 3 subcarriers
    
    pilot_syms = [6, 7, 13, 14]; % Last 2 symbols of each slot
    v_shift = mod(CellID, 6);
    
    for s = 1:N_sym
        if ismember(s, pilot_syms)
            % Base indices for pilots
            idx_p0 = [1, 7] + v_shift;
            idx_p1 = [4, 10] + v_shift; % Shifted by 3 for diversity
            
            % Wrap around
            idx_p0(idx_p0 > 12) = idx_p0(idx_p0 > 12) - 12;
            idx_p1(idx_p1 > 12) = idx_p1(idx_p1 > 12) - 12;
            
            % Pilot Value
            p_val = complex(1/sqrt(2), 1/sqrt(2));
            
            % Assign to Grids (Orthogonal Muting)
            % Antenna 0 gets Pilot at P0, Null at P1
            grid_ant0(idx_p0, s) = p_val;
            grid_ant0(idx_p1, s) = 0; 
            
            % Antenna 1 gets Pilot at P1, Null at P0
            grid_ant1(idx_p1, s) = p_val;
            grid_ant1(idx_p0, s) = 0;
            
            % Mark both as reserved so we don't put data there
            reserved_mask(idx_p0, s) = true;
            reserved_mask(idx_p1, s) = true;
        end
    end
    
    %% 3. SFBC Data Encoding (Alamouti in Frequency)
    % We process data in PAIRS.
    % Pair: [s1, s2]
    % Ant0 sends: [s1, s2]   at freq [k, k+1]
    % Ant1 sends: [-s2*, s1*] at freq [k, k+1]
    
    % Count available pairs
    available_res = (N_subc * N_sym) - sum(reserved_mask(:));
    num_pairs = floor(available_res / 2);
    
    % Generate QPSK Symbols
    M = 4; k = 2;
    total_bits_needed = num_pairs * 2 * k; % 2 symbols/pair * 2 bits/symbol
    
    if length(bits) < total_bits_needed
        bits = [bits; randi([0 1], total_bits_needed - length(bits), 1)];
    else
        bits = bits(1:total_bits_needed);
    end
    
    qam_syms = qammod(bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    % Map to Grids
    sym_idx = 1;
    
    % Iterate through grid to find available pairs of subcarriers
    for s = 1:N_sym
        for f = 1:2:11 % Step by 2 (checking f and f+1)
            % Check if f and f+1 are both free (not pilots)
            if ~reserved_mask(f, s) && ~reserved_mask(f+1, s)
                if sym_idx < length(qam_syms)
                    % Get the Pair
                    s1 = qam_syms(sym_idx);
                    s2 = qam_syms(sym_idx+1);
                    
                    % --- SFBC ENCODING ---
                    % Antenna 0
                    grid_ant0(f, s)   = s1;
                    grid_ant0(f+1, s) = s2;
                    
                    % Antenna 1 (Alamouti)
                    grid_ant1(f, s)   = -conj(s2);
                    grid_ant1(f+1, s) = conj(s1);
                    
                    sym_idx = sym_idx + 2;
                end
            end
        end
    end
end