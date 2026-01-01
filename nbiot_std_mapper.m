function [resource_grid, pilot_mask] = nbiot_std_mapper(bits, N_sym, CellID)
% NBIOT_STD_MAPPER Maps bits using Standard NB-IoT Pilot Patterns (NRS)
%
% Inputs:
%   bits    : Data bits (Column vector)
%   N_sym   : Number of symbols (Usually 14 for 1 Subframe)
%   CellID  : Physical Cell ID (Integer, determines pilot frequency shift)
%
% Outputs:
%   resource_grid : 12 x N_sym complex grid
%   pilot_mask    : Logical mask of pilot locations

    %% 1. Setup
    N_subc = 12; 
    resource_grid = zeros(N_subc, N_sym);
    pilot_mask = false(N_subc, N_sym);
    
    %% 2. Define Standard NRS Pilot Pattern
    % Rules based on 3GPP TS 36.211 (Simplified for Simulation):
    % 1. Pilots appear in the last 2 symbols of each slot.
    %    (For a 14-symbol subframe, these are indices 6, 7, 13, 14).
    % 2. Frequency separation is 6 subcarriers.
    % 3. Vertical Shift (v_shift) depends on CellID.
    
    pilot_symbols = [6, 7, 13, 14]; % 1-based indices
    v_shift = mod(CellID, 6);       % Shift can be 0 to 5
    
    for s = 1:N_sym
        if ismember(s, pilot_symbols)
            % Determine pilot subcarriers for this symbol
            % Basic pattern: indices [1, 7] shifted by v_shift
            
            % Create base indices (1-based)
            p_subs = [1, 7]; 
            
            % Apply shift
            p_subs = p_subs + v_shift;
            
            % Wrap around if shift pushes index beyond 12
            % (Note: In real LTE, it's modulo, but since BW is fixed 12, 
            % we wrap 13->1, 14->2, etc.)
            p_subs(p_subs > 12) = p_subs(p_subs > 12) - 12;
            
            % Mark grid
            pilot_mask(p_subs, s) = true;
        end
    end
    
    %% 3. Insert Pilots
    % Standard NRS pilots are complex values derived from CellID.
    % For simulation, we use a robust QPSK-like pilot: (1+j)/sqrt(2)
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
        padding = randi([0 1], max_bits - length(bits), 1);
        bits = [bits; padding];
    end
    
    % Modulate Data
    data_syms = qammod(bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    % Map to Grid (Filling only non-pilot spots)
    resource_grid(~pilot_mask) = data_syms;

end