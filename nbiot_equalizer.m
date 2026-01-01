function [rx_grid_eq, H_est_full] = nbiot_equalizer(rx_grid, CellID)
% NBIOT_EQUALIZER Performs LS Channel Estimation and Zero-Forcing Equalization
%
% Inputs:
%   rx_grid   : 12 x N_sym Received Grid (Data + Pilots + Noise)
%   CellID    : Needed to find where the pilots are hidden
%
% Outputs:
%   rx_grid_eq : Equalized Data Grid (Soft symbols)
%   H_est_full : Estimated Channel Matrix (12 x N_sym)

    [N_subc, N_sym] = size(rx_grid);
    
    %% 1. Identify Pilot Locations & Extract Values
    % We need to know EXACTLY where the pilots are to estimate H.
    % (This logic must match your nbiot_std_mapper)
    
    pilot_rows = [];
    pilot_cols = [];
    rx_pilot_vals = [];
    
    % Expected Pilot Value (Must match Transmitter)
    known_pilot = complex(1/sqrt(2), 1/sqrt(2));
    
    pilot_symbols = [6, 7, 13, 14]; % Standard NRS symbols
    v_shift = mod(CellID, 6);
    
    for s = 1:N_sym
        if ismember(s, pilot_symbols)
            % Calculate Row Indices for this symbol
            p_subs = [1, 7] + v_shift;
            p_subs(p_subs > 12) = p_subs(p_subs > 12) - 12; % Wrap
            
            % Store Coordinates
            pilot_rows = [pilot_rows; p_subs(:)];
            pilot_cols = [pilot_cols; repmat(s, length(p_subs), 1)];
            
            % Extract Received Pilot Values
            rx_pilot_vals = [rx_pilot_vals; rx_grid(p_subs, s)];
        end
    end
    
    %% 2. Least Squares (LS) Estimation
    % H_est = Y_rx / X_tx
    % Since we know X_tx is (1+j)/sqrt(2), we divide received by known.
    H_ls_points = rx_pilot_vals ./ known_pilot;
    
    %% 3. Interpolation (The "Magic" Step)
    % We have H estimates at a few scattered points. We need a full 12x14 grid.
    % We use scatteredInterpolant to fit a surface over the pilots.
    
    % Interpolator Setup
    % 'linear' fills the space between pilots.
    % 'nearest' fills the edges (e.g., symbols 1-5 which have no pilots).
    F = scatteredInterpolant(pilot_rows, pilot_cols, H_ls_points, 'linear', 'nearest');
    
    % Create a query grid for the whole frame (1 to 12, 1 to N_sym)
    [R_query, C_query] = ndgrid(1:N_subc, 1:N_sym);
    
    % Evaluate interpolator to get the Full Channel Estimate
    H_est_full = F(R_query, C_query);
    
    %% 4. Zero-Forcing Equalization
    % Divide received signal by the estimated channel
    % Y_eq = Y_rx ./ H_est
    
    % Avoid division by zero (though unlikely with fading)
    div_threshold = 1e-10; 
    H_est_safe = H_est_full;
    H_est_safe(abs(H_est_safe) < div_threshold) = div_threshold;
    
    rx_grid_eq = rx_grid ./ H_est_safe;

end