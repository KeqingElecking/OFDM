% FIXED_CHECK_RX.m
clear; clc;

%% 1. Parameters
N_fft = 64; 
N_cp = 12; 
N_sym = 14; 
CellID = 1; 
M = 4; k = log2(M); % QPSK = 2 bits

%% 2. Calculate Exact Data Capacity (To avoid Size Mismatch)
% We need to know how many REs are left for data after Pilots take their spots.
% This logic mimics the mapper's pilot allocation.
num_pilots = 0;
pilot_symbols = [6, 7, 13, 14]; % Standard NRS locations

for s = 1:N_sym
    % Logic: If this symbol has pilots, it has 2 pilots.
    % (Standard NRS uses 2 pilots per symbol in the pilot-carrying symbols)
    if ismember(mod(s-1, 14)+1, pilot_symbols) % Modulo to handle multiple subframes
        num_pilots = num_pilots + 2;
    end
end

total_REs = 12 * N_sym;
data_REs = total_REs - num_pilots;
required_bits = data_REs * k; 

fprintf('Grid Capacity: %d REs\n', total_REs);
fprintf('Pilots Used:   %d REs\n', num_pilots);
fprintf('Data Capacity: %d REs (%d bits)\n', data_REs, required_bits);

%% 3. Generate EXACT Data Bits
% Don't use gen_ofdm_data blindly, generate exactly what fits.
tx_bits = randi([0 1], required_bits, 1);

%% 4. Transmitter Chain
[grid_tx, ~] = nbiot_std_mapper(tx_bits, N_sym, CellID);
tx_wave = nbiot_ofdm_modulator(grid_tx, N_fft, N_cp);

%% 5. Channel (Optional)
USE_CHANNEL = true; 

if USE_CHANNEL
    [rx_wave, ~] = nbiot_sos_channel(tx_wave, N_fft, N_cp, 30);
else
    rx_wave = tx_wave; 
end

%% 6. Receiver Chain
[rx_bits, rx_grid] = nbiot_rx_simple(rx_wave, N_fft, N_cp, N_sym, CellID);

%% 7. Verification
% Now sizes should match exactly!
if length(tx_bits) ~= length(rx_bits)
    warning('Size mismatch! Tx: %d, Rx: %d. Truncating for BER.', length(tx_bits), length(rx_bits));
    min_len = min(length(tx_bits), length(rx_bits));
    tx_bits = tx_bits(1:min_len);
    rx_bits = rx_bits(1:min_len);
end

[num_err, ber] = biterr(tx_bits, rx_bits);

% Visualization
figure;
subplot(1,2,1);
plot(rx_grid(:), 'b.'); title('Rx Constellation');
subplot(1,2,2);
bar([0 1], [0 ber]); title(['BER: ' num2str(ber)]);

fprintf('Bit Errors: %d\n', num_err);