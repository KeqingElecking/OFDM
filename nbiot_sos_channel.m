function [rx_waveform, H_history] = nbiot_sos_channel(tx_waveform, N_fft, N_cp, snr_db)
% NBIOT_SOS_CHANNEL Rayleigh Multipath Channel using Sum-of-Sinusoids (Jakes)
% 
% Based on user-provided logic:
%   - N1=9, N2=10 sinusoids
%   - Deterministic Jakes frequencies
%   - Random phases
%   - Multipath (6 Taps)
%
% Inputs:
%   tx_waveform : Serial Time Domain Signal
%   N_fft, N_cp : System Parameters
%   snr_db      : Noise Level
%
% Outputs:
%   rx_waveform : Faded Signal + Noise
%   H_history   : Channel Response for debugging

    %% 1. SoS Parameters (From your code)
    f_m = 91;        % Max Doppler (Hz)
    b = 0.5;         % Variance parameter
    N1 = 9;          % Sinusoids for In-phase
    N2 = 10;         % Sinusoids for Quadrature
    
    % Power Delay Profile (6 Taps Exponential)
    L_taps = 6;
    pdp = exp(-(0:L_taps-1));
    pdp = pdp / sum(pdp); % Normalize total power to 1
    
    %% 2. Time Vector Setup
    % We need 't' to match the length of the input signal
    % Sampling Rate for NB-IoT (64 FFT) is typically 960 kHz (15kHz * 64)
    fs = 15e3 * N_fft; 
    total_samples = length(tx_waveform);
    t = (0:total_samples-1) / fs;
    
    %% 3. Generate Fading Processes (One per Tap)
    % We generate 'L_taps' independent fading waveforms
    h_taps_time = zeros(L_taps, total_samples);
    
    for tap = 1:L_taps
        % --- A. Jakes' Setup (Your Logic) ---
        c1 = zeros(1, N1); f1 = zeros(1, N1); th1 = zeros(1, N1);
        c2 = zeros(1, N2); f2 = zeros(1, N2); th2 = zeros(1, N2);
        
        % In-Phase Constants
        for n = 1:N1
            c1(n) = sqrt(2*b/N1);
            f1(n) = f_m * sin(pi*(n-0.5)/(2*N1));
            th1(n) = rand(1) * 2 * pi; % Random Phase
        end
        
        % Quadrature Constants
        for n = 1:N2
            c2(n) = sqrt(2*b/N2);
            f2(n) = f_m * sin(pi*(n-0.5)/(2*N2));
            th2(n) = rand(1) * 2 * pi; % Random Phase
        end
        
        % --- B. Sum of Sinusoids (Vectorized 'g' function) ---
        % g1 = sum( c1 * cos(2*pi*f1*t + th1) )
        g1 = zeros(1, total_samples);
        for n = 1:N1
            g1 = g1 + c1(n) * cos(2*pi*f1(n)*t + th1(n));
        end
        
        g2 = zeros(1, total_samples);
        for n = 1:N2
            g2 = g2 + c2(n) * cos(2*pi*f2(n)*t + th2(n));
        end
        
        % Combine I + jQ and scale by PDP power
        % h = (g1 + 1j*g2) * sqrt(Tap_Power)
        h_taps_time(tap, :) = (g1 + 1j*g2) * sqrt(pdp(tap));
    end
    
    %% 4. Apply Multipath Convolution (Time-Varying)
    % Since h changes every sample, we use Overlap-Add logic
    
    rx_waveform = zeros(total_samples + L_taps, 1);
    
    % To speed up, we assume block fading per Symbol (approximate)
    % or sample-by-sample convolution. 
    % Here is a precise sample-by-sample filter implementation:
    
    padded_tx = [zeros(L_taps-1, 1); tx_waveform];
    
    for n = 1:total_samples
        % The received sample y[n] is the sum of delayed inputs x[n-k] * h_k[n]
        % y[n] = sum( x[n-k] * h[k, n] )
        
        for k = 1:L_taps
            if (n-k+1 > 0)
                val_tx = tx_waveform(n-k+1);
                val_h  = h_taps_time(k, n); % Instantaneous tap value
                rx_waveform(n) = rx_waveform(n) + val_tx * val_h;
            end
        end
    end
    
    rx_waveform = rx_waveform(1:total_samples);
    
    %% 5. Output Data
    % Save H_history (sampled once per symbol for the Equalizer)
    sym_len = N_fft + N_cp;
    num_syms = floor(total_samples / sym_len);
    H_history = zeros(L_taps, num_syms);
    
    % Downsample H to symbol rate (taking the middle of the symbol)
    for s = 1:num_syms
        idx = round((s-1)*sym_len + sym_len/2);
        H_history(:, s) = h_taps_time(:, idx);
    end
    
    % Add Noise
    rx_waveform = awgn(rx_waveform, snr_db, 'measured');

end