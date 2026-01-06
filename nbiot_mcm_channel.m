function [rx_waveform, H_history] = nbiot_mcm_channel(tx_waveform, N_fft, N_cp, snr_db)
% NBIOT_MCM_CHANNEL Applies Time-Varying Multipath Fading (Sum-of-Sinusoids)
%
% Inputs:
%   tx_waveform : Serial time-domain signal (from nbiot_ofdm_modulator)
%   N_fft       : FFT Size (e.g., 64)
%   N_cp        : Cyclic Prefix Length (e.g., 16)
%   snr_db      : SNR in dB for AWGN
%
% Outputs:
%   rx_waveform : Received serial signal (with Fading + ISI + Noise)
%   H_history   : Matrix [N_taps x N_sym] of the channel Impulse Response 
%                 (Used for Perfect Channel Estimation / Debugging)

    %% 1. Channel Physics Parameters (Derived from your legacy code)
    % Power Delay Profile (approximating 'rho.am')
    % We use a standard Exponential Decay profile for multipath
    L_taps = 6;
    rho = exp(-(0:L_taps-1)); 
    rho = rho / norm(rho); % Normalize energy to 1
    
    % Sum-of-Sinusoids Parameters
    num_summations = 40; % As per your code
    f_dmax = 150.0;       % Max Doppler Shift (Hz)
    
    % Time Parameters
    % Calculate symbol duration based on NB-IoT 15kHz spacing
    % T_sym = 1/15kHz (approx 66.7us) + CP. 
    % For simulation, we just need a step size 't' relative to f_dmax.
    sym_duration = (N_fft + N_cp) * 50e-9; % Using your legacy t_a = 50ns
    
    %% 2. Initialization
    % Calculate number of symbols in the input stream
    sym_len = N_fft + N_cp;
    num_syms = length(tx_waveform) / sym_len;
    
    if floor(num_syms) ~= num_syms
        error('Input waveform length must be a multiple of (N_fft + N_cp)');
    end
    
    % Initialize Random Phases (Fixed for the whole duration to maintain correlation)
    % This is 'u' in your legacy code.
    u = rand(L_taps, num_summations);
    
    % Prepare Output Buffer (Length + tail for last symbol ISI)
    out_len = length(tx_waveform) + L_taps - 1;
    rx_stream_accum = zeros(out_len, 1);
    
    % History buffer for channel estimation
    H_history = zeros(L_taps, num_syms);
    
    t_current = 0;
    
    %% 3. Processing Loop (Symbol-by-Symbol)
    for i = 1:num_syms
        % A. Extract current OFDM symbol (Time Domain)
        idx_start = (i-1)*sym_len + 1;
        idx_end   = i*sym_len;
        sym_tx = tx_waveform(idx_start:idx_end);
        
        % B. Generate Channel Impulse Response 'h' for this moment
        % (MCM_channel_model logic from your legacy code)
        h_inst = zeros(L_taps, 1);
        
        for k = 1:L_taps
            u_k = u(k, :);           % Random variable for kth tap
            phi = 2 * pi * u_k;      % Phase coefficients
            
            % Doppler frequency for each scatterer
            f_d = f_dmax * sin(2 * pi * u_k); 
            
            % Sum-of-Sinusoids Formula
            % Summing complex exponentials to get Rayleigh fading coefficient
            scatterers = exp(1j * phi) .* exp(1j * 2 * pi * f_d * t_current);
            h_inst(k) = rho(k) * (1/sqrt(num_summations)) * sum(scatterers);
        end
        
        % Store h for debugging/receiver
        H_history(:, i) = h_inst;
        
        % C. Apply Multipath (Convolution)
        % This creates the signal + tail (ISI)
        sym_rx_raw = conv(sym_tx, h_inst);
        
        % D. Overlap-Add to main stream
        % We add this symbol's response to the stream. 
        % The "tail" of this symbol will overlap with the "head" of the next.
        idx_out_end = idx_start + length(sym_rx_raw) - 1;
        rx_stream_accum(idx_start:idx_out_end) = rx_stream_accum(idx_start:idx_out_end) + sym_rx_raw;
        
        % E. Update Time
        t_current = t_current + sym_duration;
    end
    
    % Truncate to original length (optional, usually Rx captures tails)
    rx_waveform_clean = rx_stream_accum(1:length(tx_waveform));
    
    %% 4. Add Noise (AWGN)
    rx_waveform = awgn(rx_waveform_clean, snr_db, 'measured');

end