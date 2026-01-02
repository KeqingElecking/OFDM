function [out_bits] = nbiot_channel_coding(in_bits, type, direction)
% NBIOT_CHANNEL_CODING Implements simple FEC schemes
% Inputs:
%   in_bits   : Input bit vector
%   type      : 'Hamming74' or 'Repetition3'
%   direction : 'ENCODE' or 'DECODE'

    switch type
        case 'Hamming74'
            % --- DEFINE MATRIX G GLOBALLY FOR THIS CASE ---
            % We need G for both Encoding and Decoding (valid code search)
            G = [1 1 0 1; 
                 1 0 1 1; 
                 1 0 0 0; 
                 0 1 1 1; 
                 0 1 0 0; 
                 0 0 1 0; 
                 0 0 0 1];

            if strcmp(direction, 'ENCODE')
                % Pad input to multiple of 4
                n = length(in_bits);
                rem_bits = mod(n, 4);
                if rem_bits > 0, in_bits = [in_bits; zeros(4-rem_bits, 1)]; end
                
                % Reshape to 4-bit chunks
                msg = reshape(in_bits, 4, []).';
                
                % Encode: c = m * G' (modulo 2)
                code = mod(msg * G', 2);
                out_bits = reshape(code.', [], 1);
                
            else % DECODE
                % Reshape to 7-bit chunks
                n = length(in_bits);
                % Truncate to multiple of 7 (safety against padding errors)
                n_clean = floor(n/7)*7;
                code = reshape(in_bits(1:n_clean), 7, []).';
                
                % --- MINIMUM DISTANCE DECODER ---
                % Generate all 16 valid codewords using G
                valid_msgs = dec2bin(0:15, 4) - '0';
                
                % Now G is accessible here!
                valid_codes = mod(valid_msgs * G', 2);
                
                decoded_msg = zeros(size(code,1), 4);
                
                for i = 1:size(code, 1)
                    rx_word = code(i, :);
                    
                    % Calculate Hamming distance to all 16 valid codes
                    % dist = sum( abs(Valid - Rx) )
                    dist = sum(abs(valid_codes - rx_word), 2);
                    
                    % Find the closest match
                    [~, min_idx] = min(dist);
                    decoded_msg(i, :) = valid_msgs(min_idx, :);
                end
                
                out_bits = reshape(decoded_msg.', [], 1);
            end

        case 'Repetition3'
            if strcmp(direction, 'ENCODE')
                % Repeat every bit 3 times: 1 -> 1 1 1
                out_bits = repelem(in_bits, 3);
            else
                % Decode: Majority Vote
                % Reshape to columns of 3
                n = length(in_bits);
                n_clean = floor(n/3)*3;
                mat = reshape(in_bits(1:n_clean), 3, []);
                
                % Sum columns. If sum >= 2, bit is 1. Else 0.
                votes = sum(mat, 1);
                out_bits = (votes >= 2).'; 
            end
            
        otherwise
            out_bits = in_bits; % Pass through
    end
end