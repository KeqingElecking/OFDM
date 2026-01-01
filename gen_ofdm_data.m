function [tx_bits, tx_syms] = gen_ofdm_data(N_data_subc, N_sym, M, N_tx)
% GEN_OFDM_DATA Generates random data for a MIMO-OFDM simulation
%
% Inputs:
%   N_data_subc : Number of active data subcarriers (e.g., 52 or 48)
%   N_sym       : Number of OFDM time symbols to simulate
%   M           : Modulation Order (e.g., 4=QPSK, 16=16QAM, 64=64QAM)
%   N_tx        : Number of Transmit Antennas (for MIMO Spatial Mux)
%                 *Set to 1 for SISO*
%
% Outputs:
%   tx_bits     : Column vector of binary data (0s and 1s)
%   tx_syms     : Complex QAM symbols reshaped for the grid
%                 Dimensions: [N_data_subc, N_sym, N_tx]

    % 1. Calculate parameters
    k_bits = log2(M); % Bits per symbol (e.g., 64QAM -> 6 bits)
    
    % Total number of QAM symbols needed across all time, freq, and antennas
    total_qam_syms = N_data_subc * N_sym * N_tx;
    
    % Total number of bits needed
    total_bits = total_qam_syms * k_bits;
    
    % 2. Generate Random Bits
    % We use randi([0 1]) to simulate a binary data stream
    tx_bits = randi([0 1], total_bits, 1);
    
    % 3. Modulate (Bit -> Complex Symbol)
    % UnitAveragePower=true ensures the constellation has energy = 1
    % InputType='bit' tells MATLAB to group the bits automatically
    mod_data = qammod(tx_bits, M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    % 4. Reshape for the OFDM Grid
    % We organize this into a 3D Matrix: [Frequency x Time x Antenna]
    tx_syms = reshape(mod_data, N_data_subc, N_sym, N_tx);

end