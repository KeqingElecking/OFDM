function [grid_ant0, grid_ant1] = nbiot_dl_spatial_mapper(bits, N_sym, CellID)
% Maps bits to 2 Layers (Spatial Multiplexing) for Downlink OFDM
    
    N_subc = 12;
    grid_ant0 = zeros(N_subc, N_sym);
    grid_ant1 = zeros(N_subc, N_sym);
    
    %% 1. Pilot Setup (Ports 0 & 1 Orthogonal)
    pilot_mask = false(N_subc, N_sym);
    pilot_syms = [];
    
    % Standard NRS locations: Symbols 5, 6 (0-based) -> 6, 7 (1-based)
    for s = 1:N_sym
        if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
            pilot_syms = [pilot_syms, s];
        end
    end
    
    v_shift = mod(CellID, 6);
    p_val = complex(1/sqrt(2), 1/sqrt(2));
    
    for s = 1:N_sym
        if ismember(s, pilot_syms)
            % Antenna 0 Indices
            idx_p0 = [1, 7] + v_shift;
            idx_p0(idx_p0 > 12) = idx_p0(idx_p0 > 12) - 12;
            
            % Antenna 1 Indices
            idx_p1 = [4, 10] + v_shift;
            idx_p1(idx_p1 > 12) = idx_p1(idx_p1 > 12) - 12;
            
            % Ant 0 transmits P0, mutes P1
            grid_ant0(idx_p0, s) = p_val;
            grid_ant0(idx_p1, s) = 0;
            
            % Ant 1 transmits P1, mutes P0
            grid_ant1(idx_p1, s) = p_val;
            grid_ant1(idx_p0, s) = 0;
            
            pilot_mask(idx_p0, s) = true;
            pilot_mask(idx_p1, s) = true;
        end
    end
    
    %% 2. Data Splitting & Mapping
    % Calculate Capacity
    num_data_res = sum(~pilot_mask(:));
    M = 4; k = 2; % QPSK
    max_bits = num_data_res * k * 2; % x2 for two layers
    
    % Pad/Truncate
    if length(bits) < max_bits
        bits = [bits; randi([0 1], max_bits - length(bits), 1)];
    else
        bits = bits(1:max_bits);
    end
    
    % Split
    bits_L1 = bits(1:max_bits/2);
    bits_L2 = bits(max_bits/2+1:end);
    
    % Map
    syms_L1 = qammod(bits_L1, M, 'InputType', 'bit', 'UnitAveragePower', true);
    syms_L2 = qammod(bits_L2, M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    % Fill Grid
    grid_ant0(~pilot_mask) = syms_L1;
    grid_ant1(~pilot_mask) = syms_L2;
end