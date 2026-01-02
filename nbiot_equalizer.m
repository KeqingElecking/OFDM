function [rx_grid_eq, H_est_full] = nbiot_equalizer(rx_grid, CellID)
% NBIOT_EQUALIZER Performs LS Channel Estimation and Zero-Forcing Equalization
%
% FIX APPLIED: Dynamic Pilot Detection for arbitrary N_sym length

    [N_subc, N_sym] = size(rx_grid);
    
    %% 1. Identify Pilot Locations & Extract Values (FIXED)
    pilot_rows = [];
    pilot_cols = [];
    rx_pilot_vals = [];
    
    % Expected Pilot Value (Must match Transmitter)
    known_pilot = complex(1/sqrt(2), 1/sqrt(2));
    
    v_shift = mod(CellID, 6);
    
    % Loop through ALL symbols to find pilots dynamically
    for s = 1:N_sym
        % Standard NB-IoT: Pilots at last 2 symbols of every slot (indices 5,6)
        if mod(s-1, 7) == 5 || mod(s-1, 7) == 6
            
            % Calculate Row Indices for this symbol
            p_subs = [1, 7] + v_shift;
            p_subs(p_subs > 12) = p_subs(p_subs > 12) - 12; % Wrap
            
            % Store Coordinates for Interpolation
            pilot_rows = [pilot_rows; p_subs(:)];
            pilot_cols = [pilot_cols; repmat(s, length(p_subs), 1)];
            
            % Extract Received Pilot Values
            rx_pilot_vals = [rx_pilot_vals; rx_grid(p_subs, s)];
        end
    end
    
    %% 2. Least Squares (LS) Estimation
    % H_est = Y_rx / X_tx
    H_ls_points = rx_pilot_vals ./ known_pilot;
    
    %% 3. Interpolation
    % Use scatteredInterpolant to fit a surface over the pilots.
    if isempty(pilot_rows)
        % Fallback if no pilots found (e.g., extremely short sim)
        H_est_full = ones(N_subc, N_sym);
    else
        F = scatteredInterpolant(pilot_rows, pilot_cols, H_ls_points, 'linear', 'nearest');
        
        % Create a query grid for the whole frame
        [R_query, C_query] = ndgrid(1:N_subc, 1:N_sym);
        
        % Evaluate interpolator
        H_est_full = F(R_query, C_query);
    end
    
    %% 4. Zero-Forcing Equalization
    % Y_eq = Y_rx ./ H_est
    div_threshold = 1e-10; 
    H_est_safe = H_est_full;
    H_est_safe(abs(H_est_safe) < div_threshold) = div_threshold;
    
    rx_grid_eq = rx_grid ./ H_est_safe;

end