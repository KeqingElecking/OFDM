function [grid_ant0, grid_ant1] = nbiot_ul_spatial_mapper(bits, N_sym, M)
% NBIOT_UL_SPATIAL_MAPPER Maps bits to 2 Layers for Uplink SC-FDMA
% FIX: Pre-codes pilots with IFFT so they are correct in Frequency Domain
    
    N_subc = 12;
    grid_ant0 = zeros(N_subc, N_sym);
    grid_ant1 = zeros(N_subc, N_sym);
    pilot_mask = false(N_subc, N_sym);
    
    %% 1. Pilot Setup (Symbol 3)
    p_val = complex(1/sqrt(2), 1/sqrt(2));
    
    % Define Frequency Domain Pilot Combs
    % Ant 0: Active on 1, 3, 5... (Odd)
    P_freq_0 = zeros(N_subc, 1);
    P_freq_0(1:2:12) = p_val;
    
    % Ant 1: Active on 2, 4, 6... (Even)
    P_freq_1 = zeros(N_subc, 1);
    P_freq_1(2:2:12) = p_val;
    
    % --- THE FIX: Pre-code to Time Domain ---
    % Because SC-FDMA Modulator does FFT, we must do IFFT here.
    p_time_0 = ifft(P_freq_0);
    p_time_1 = ifft(P_freq_1);
    % ----------------------------------------
    
    for s = 1:N_sym
        if mod(s-1, 7) == 3 
            % Place Time-Domain sequences into the grid
            grid_ant0(:, s) = p_time_0;
            grid_ant1(:, s) = p_time_1;
            
            % Mark specific subcarriers as pilots for the bit loading count
            % (Conceptually, all are occupied by the sequence)
            pilot_mask(:, s) = true; 
        end
    end
    
    %% 2. Data Splitting
    num_data_res = sum(~pilot_mask(:));
    k = log2(M);
    max_bits = num_data_res * k * 2; 
    
    if length(bits) < max_bits
        bits = [bits; randi([0 1], max_bits - length(bits), 1)];
    else
        bits = bits(1:max_bits);
    end
    
    bits_L1 = bits(1:max_bits/2);
    bits_L2 = bits(max_bits/2+1:end);
    
    syms_L1 = qammod(bits_L1, M, 'InputType', 'bit', 'UnitAveragePower', true);
    syms_L2 = qammod(bits_L2, M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    grid_ant0(~pilot_mask) = syms_L1;
    grid_ant1(~pilot_mask) = syms_L2;
end