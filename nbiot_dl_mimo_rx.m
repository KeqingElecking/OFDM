function [rx_bits, grid_eq] = nbiot_dl_mimo_rx(rx_wave, N_fft, N_cp, N_sym, CellID)
% NBIOT_DL_MIMO_RX Downlink Receiver for SFBC (2 Tx -> 1 Rx)

    N_subc = 12;
    
    %% 1. OFDM Demodulation
    [~, rx_grid] = nbiot_rx_simple(rx_wave, N_fft, N_cp, N_sym, CellID);
    
    %% 2. MIMO Channel Estimation
    % Dynamic Pilot Logic (Matches nbiot_std_mapper)
    pilot_symbols = [];
    for n = 1:N_sym
        % Check 6th and 7th symbol of every slot
        if mod(n-1, 7) == 5 || mod(n-1, 7) == 6
            pilot_symbols = [pilot_symbols, n];
        end
    end 
    
    v_shift = mod(CellID, 6);
    mask_p1 = false(12, N_sym);
    mask_p2 = false(12, N_sym);
    
    for s = 1:N_sym
        % FIX 1: Use the calculated 'pilot_symbols' variable
        if ismember(s, pilot_symbols)
            % Ant 1 Indices
            idx1 = [1, 7] + v_shift;
            idx1(idx1>12) = idx1(idx1>12)-12;
            mask_p1(idx1, s) = true;
            
            % Ant 2 Indices (Shift +3)
            idx2 = [4, 10] + v_shift;
            idx2(idx2>12) = idx2(idx2>12)-12;
            mask_p2(idx2, s) = true;
        end
    end
    
    % Estimator Helper
    known_p = complex(1/sqrt(2), 1/sqrt(2));
    H1_est = interp_channel(rx_grid, mask_p1, known_p);
    H2_est = interp_channel(rx_grid, mask_p2, known_p);
    
    %% 3. SFBC Decoding
    grid_eq = zeros(12, N_sym);
    all_pilots = mask_p1 | mask_p2;
    data_mask = false(12, N_sym); % FIX 3: Track valid data locations
    
    for s = 1:N_sym
        for f = 1:2:11
            % FIX 2: Check BOTH subcarriers in the pair. 
            % If either is a pilot, SFBC cannot be performed for this pair.
            if ~all_pilots(f, s) && ~all_pilots(f+1, s)
                
                r1 = rx_grid(f, s);
                r2 = rx_grid(f+1, s);
                
                h1 = H1_est(f, s);
                h2 = H2_est(f, s);
                
                % Alamouti Combiner
                s1_hat = conj(h1)*r1 + h2*conj(r2);
                s2_hat = conj(h2)*r1 - h1*conj(r2);
                
                pow_norm = abs(h1)^2 + abs(h2)^2 + 1e-10;
                grid_eq(f, s)   = s1_hat / pow_norm;
                grid_eq(f+1, s) = s2_hat / pow_norm;
                
                % Mark these REs as containing valid data
                data_mask(f, s) = true;
                data_mask(f+1, s) = true;
            end
        end
    end
    
    %% 4. Demodulate
    % Only extract REs where SFBC was actually performed
    rx_data = grid_eq(data_mask);
    rx_bits = qamdemod(rx_data, 4, 'OutputType', 'bit', 'UnitAveragePower', true);
end

function H_full = interp_channel(grid, mask, p_val)
    [r, c] = find(mask);
    if isempty(r), H_full = ones(size(grid)); return; end
    h_vals = grid(mask) ./ p_val;
    F = scatteredInterpolant(r, c, h_vals, 'linear', 'nearest');
    [R, C] = ndgrid(1:12, 1:size(grid,2));
    H_full = F(R, C);
end