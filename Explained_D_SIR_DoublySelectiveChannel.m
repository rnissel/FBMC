% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script calculates the Signal-to-Interference Ratio (SIR)
% It includes different methods.
%       1) The matrix description, see Section IV.A and IV.B. 
%       2) Based on our modulation objects, which is similar to the matrix
%          description but more efficient.
%       3) We include a Monte Carlo simulation of the SIR in order to check
%          the theory
%       4) The theoretical solution based on the generalized hypergeometric
%          function. Valid for CP-OFDM, a Jakes Doppler spectrum, no ISI
%          and infinitely many subcarriers.
% All this SIR calulations/simulations should lead to approximatly the same
% result. The matrix calculation is very inefficient => we choose a small
% sampling rate and only a few time-frequency positions.

clear; close all;

%% Parameters
Channel_Velocity_kmh            = 250;                                  % Velocity of the user
Channel_CarrierFrequency        = 2.5e9;                                % Carrier frequency
Channel_PowerDelayProfile       = [1 0.1 0.01];                         % Represents the Power delay profile. We define it sample wise and not by, e.g., "PedestrianA" because the sampling rate is too low to accuratly describe practical channel models. The main purpose here is to show how the calculation works in principle, whereas we want to keep the computational complexity as low as possible.

Simulation_MonteCarloRepetition = 10000;

%% Modulation objects
% We choose fixed values here just as an example.
% FBMC Object
FBMC = Modulation.FBMC(...
    12,...                                                              % Number of subcarriers
    8,...                                                               % Number FBMC symbols in time
    15e3,...                                                            % Subcarrier spacing (Hz)
    15e3*12,...                                                         % Sampling rate (Samples/s), critically sampled to investigate the effect of infinitely many subcarriers
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz).  Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    'Hermite-OQAM',...                                                  % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM. The data rate of QAM is reduced by a factor of two compared to OQAM, but robustness in doubly-selective channels is inceased
    8, ...                                                              % Overlapping factor (also determines oversampling in the frequency domain)                                   
    0, ...                                                              % Initial phase shift
    true ...                                                            % Efficient IFFT implementation
    );
FBMC_BlockOverlapTime = (FBMC.PrototypeFilter.OverlappingFactor-1/2)*FBMC.PHY.TimeSpacing;

% OFDM Object (without CP)
OFDMnoCP = Modulation.OFDM(...
    FBMC.Nr.Subcarriers,...                                             % Number of subcarriers
    FBMC.Nr.MCSymbols/2,...                                             % Number OFDM symbols in time                                                 
    FBMC.PHY.SubcarrierSpacing,...                                      % Subcarrier spacing (Hz) 
    FBMC.PHY.SamplingRate,...                                           % Sampling rate (Samples/s), critically sampled to investigate the effect of infinitly many subcarriers                                       
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    0, ...                                                              % Length of the cyclic prefix (s)                 
    FBMC_BlockOverlapTime ...                                           % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    );

% CP-OFDM Object
OFDM = Modulation.OFDM(...
    FBMC.Nr.Subcarriers,...                                             % Number of subcarriers
    FBMC.Nr.MCSymbols/4,...                                             % Number OFDM symbols in time                                                 
    FBMC.PHY.SubcarrierSpacing,...                                      % Subcarrier spacing (Hz) 
    FBMC.PHY.SamplingRate,...                                           % Sampling rate (Samples/s), critically sampled to investigate the effect of infinitly many subcarriers                                       
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    1/FBMC.PHY.SubcarrierSpacing, ...                                   % Length of the cyclic prefix (s)                
    FBMC_BlockOverlapTime ...                                           % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    );

% Channel model object
ChannelModel = Channel.FastFading(...
    FBMC.PHY.SamplingRate,...                                           % Sampling rate (Samples/s)
    Channel_PowerDelayProfile,...                                       % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
    FBMC.Nr.SamplesTotal,...                                            % Number of total samples
    Channel_Velocity_kmh/3.6*Channel_CarrierFrequency/2.998e8,...       % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    'Jakes',...                                                         % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
    200, ...                                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
    1,...                                                               % Number of transmit antennas
    1,...                                                               % Number of receive antennas
    true ...                                                            % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );

%% Calculate transmit and receive matrices for FBMC and OFDM
G_TX_FBMC       = FBMC.GetTXMatrix;
G_RX_FBMC       = FBMC.GetRXMatrix';
G_TX_OFDM       = OFDM.GetTXMatrix;
G_RX_OFDM       = OFDM.GetRXMatrix';
G_TX_OFDMnoCP   = OFDMnoCP.GetTXMatrix;
G_RX_OFDMnoCP   = OFDMnoCP.GetRXMatrix';

%% Vectorized channel correlation matrix
R_vecH = ChannelModel.GetCorrelationMatrix;


%% SIR calculation in Matrix Notation, see Section IV.
% Find the middle position of the time-frequency grid in vectorized form.
% We calculate the SIR only at the middle position
PosFBMC             = round(FBMC.Nr.Subcarriers/2)      +   FBMC.Nr.Subcarriers*round(FBMC.Nr.MCSymbols/2-1);
PosOFDM             = round(OFDM.Nr.Subcarriers/2)      +   OFDM.Nr.Subcarriers*round(OFDM.Nr.MCSymbols/2-1);
PosOFDMnoCP         = round(OFDMnoCP.Nr.Subcarriers/2)  +   OFDMnoCP.Nr.Subcarriers*round(OFDMnoCP.Nr.MCSymbols/2-1);

% QAM calculation Equation (36)
Gamma_OFDM          = kron(G_TX_OFDM.',G_RX_OFDM(:,PosOFDM)')*R_vecH*kron(G_TX_OFDM.',G_RX_OFDM(:,PosOFDM)')';
Gamma_OFDMnoCP      = kron(G_TX_OFDMnoCP.',G_RX_OFDMnoCP(:,PosOFDMnoCP)')*R_vecH*kron(G_TX_OFDMnoCP.',G_RX_OFDMnoCP(:,PosOFDMnoCP)')';
% OQAM calculation
Gamma_FBMC          = kron(G_TX_FBMC.',G_RX_FBMC(:,PosFBMC)')*R_vecH*kron(G_TX_FBMC.',G_RX_FBMC(:,PosFBMC)')';
% Equation (37)
[V,D]               = eig(Gamma_FBMC);
Omega               = V*sqrt(D);      
% Equation (38)
OmegaTilde = nan(size(Omega));
for i_column = 1:size(OmegaTilde,2)
   for i_row = 1:size(OmegaTilde,1)
       OmegaTilde(i_row,i_column) = Omega(i_row,i_column)*abs(Omega(PosFBMC,i_column))/Omega(PosFBMC,i_column);
   end    
end
% Equation (39)
GammaTilde_FBMC = real(OmegaTilde)*real(OmegaTilde)';

% Calculat the SIR (35) and (40)
SIR_OFDM_Matrix     = 10*log10(abs(Gamma_OFDM(PosOFDM,PosOFDM)/(trace(Gamma_OFDM)-Gamma_OFDM(PosOFDM,PosOFDM)))); % abs due to num inaccuracies
SIR_OFDMnoCP_Matrix = 10*log10(abs(Gamma_OFDMnoCP(PosOFDMnoCP,PosOFDMnoCP)/(trace(Gamma_OFDMnoCP)-Gamma_OFDMnoCP(PosOFDMnoCP,PosOFDMnoCP)))); % abs due to num inaccuracies
SIR_FBMC_Matrix     = 10*log10(abs(GammaTilde_FBMC(PosFBMC,PosFBMC)/(trace(GammaTilde_FBMC)-GammaTilde_FBMC(PosFBMC,PosFBMC)))); % abs due to num inaccuracies


%% Built-in SIR calculation of the modulation objects (similar to the matrix representation but more efficient!)
[PS_OFDM,PI_OFDM] = OFDM.GetSignalAndInterferencePowerQAM(...
                R_vecH,...                                                  % Let the received signal be r=H*s with H representing a time-variant convolution matrix. Then "VectorizedChannelCorrelationMatrix" represents the expectation E{{H(:)*H(:)'}. We can obtain such matrix by ChannelModel.GetCorrelationMatrix
                eye(OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols),...              % Correlation matrix of the vectorized data symbols
                0,...                                                       % Time offset in samples (to improve the SIR)
                round(OFDM.Nr.Subcarriers/2),...                            % Subcarrier position for which the SIR is calculated.
                round(OFDM.Nr.MCSymbols/2)...                               % OFDM symbol position in time for which the SIR is calculated.
                );
[PS_OFDMnoCP,PI_OFDMnoCP] = OFDMnoCP.GetSignalAndInterferencePowerQAM(...
                R_vecH,...                                                  % Let the received signal be r=H*s with H representing a time-variant convolution matrix. Then "VectorizedChannelCorrelationMatrix" represents the expectation E{{H(:)*H(:)'}. We can obtain such matrix by ChannelModel.GetCorrelationMatrix
                eye(OFDMnoCP.Nr.Subcarriers*OFDMnoCP.Nr.MCSymbols),...      % Correlation matrix of the vectorized data symbols
                0,...                                                       % Time offset in samples (to improve the SIR)
                round(OFDMnoCP.Nr.Subcarriers/2),...                        % Subcarrier position for which the SIR is calculated.
                round(OFDMnoCP.Nr.MCSymbols/2)...                           % OFDM symbol position in time for which the SIR is calculated.
                );  
[PS_FBMC,PI_FBMC] = FBMC.GetSignalAndInterferencePowerOQAM(...
                R_vecH,...                                                  % Let the received signal be r=H*s with H representing a time-variant convolution matrix. Then "VectorizedChannelCorrelationMatrix" represents the expectation E{{H(:)*H(:)'}. We can obtain such matrix by ChannelModel.GetCorrelationMatrix
                eye(FBMC.Nr.Subcarriers*FBMC.Nr.MCSymbols),...              % Correlation matrix of the vectorized data symbols
                0,...                                                       % Time offset in samples (to improve the SIR)
                round(FBMC.Nr.Subcarriers/2),...                            % Subcarrier position for which the SIR is calculated.
                round(FBMC.Nr.MCSymbols/2)...                               % FBMC symbol position in time for which the SIR is calculated.
                );                          
SIR_OFDM_Object         = 10*log10( PS_OFDM      /  PI_OFDM);
SIR_OFDMnoCP_Object     = 10*log10( PS_OFDMnoCP  /  PI_OFDMnoCP);
SIR_FBMC_Object         = 10*log10( PS_FBMC      /  PI_FBMC);
    

%% Theory OFDM, Jakes Doppler spectrum, no ISI, similar to (42) but we do not consider the optimal subcarrier spacing here
% Note that this theory assumes continuous time => the SIR of our matrix 
% formulation is only valid in the limit case of sampling rate => inf. 
% However the error is so small that it can be neglected. 
PS_OFDM_Jakes           = hypergeom([1/2],[3/2 2],-(ChannelModel.PHY.MaximumDopplerShift/OFDM.PHY.SubcarrierSpacing*pi).^2);
PI_OFDM_Jakes           = 1 - PS_OFDM_Jakes;
SIR_OFDM_TheoryJakes    = 10*log10( PS_OFDM_Jakes / PI_OFDM_Jakes ) ;


%% Simulation (to check theory)
for i_rep = 1:Simulation_MonteCarloRepetition
    % Generate new channel realization
    ChannelModel.NewRealization;
    % Time-variant convolution matrix
    H = ChannelModel.GetConvolutionMatrix{1};

    % Transmission vector, i.e., how symbol y_{l,k} is affected by vector x, see (34)
    d_OFDM          = G_RX_OFDM(:,PosOFDM)'         * H * G_TX_OFDM;
    d_OFDMnoCP      = G_RX_OFDMnoCP(:,PosOFDMnoCP)' * H * G_TX_OFDMnoCP;
    d_FBMC          = G_RX_FBMC(:,PosFBMC)'         * H * G_TX_FBMC;
    % Equalize the phase in FBMC and take the real part
    dTilde_FBMC     = real(d_FBMC * abs(d_FBMC(PosFBMC)) / d_FBMC(PosFBMC) );
    
    % Signal power PS and interference power PI
    PS_OFDM_Simulation(i_rep)       = abs(d_OFDM(PosOFDM)).^2;
    PI_OFDM_Simulation(i_rep)       = sum(abs(d_OFDM(:)).^2) - PS_OFDM_Simulation(i_rep);
    
    PS_OFDMnoCP_Simulation(i_rep)   = abs(d_OFDMnoCP(PosOFDMnoCP)).^2;
    PI_OFDMnoCP_Simulation(i_rep)   = sum(abs(d_OFDMnoCP(:)).^2) - PS_OFDMnoCP_Simulation(i_rep);   
    
    PS_FBMC_Simulation(i_rep)       = abs(dTilde_FBMC(PosFBMC)).^2;
    PI_FBMC_Simulation(i_rep)       = sum(abs(dTilde_FBMC(:)).^2) - PS_FBMC_Simulation(i_rep);
   
    if mod(i_rep,100)==0
        disp(['Simulation: ' int2str(i_rep/Simulation_MonteCarloRepetition*100) '%']);
    end
end
SIR_OFDM_Simulation     = 10*log10(sum(PS_OFDM_Simulation)/sum(PI_OFDM_Simulation));
SIR_OFDMnoCP_Simulation = 10*log10(sum(PS_OFDMnoCP_Simulation)/sum(PI_OFDMnoCP_Simulation));
SIR_FBMC_Simulation     = 10*log10(sum(PS_FBMC_Simulation)/sum(PI_FBMC_Simulation));


%% Plot Results
fprintf('==============================================================================================\n');
fprintf('==== SIR comparision (matrix, object, sim. and  hyper. fun. should be (approx.) the same) ====\n');
fprintf('               | SIR (Matrix) | SIR (Object) | SIR(Simulation)| SIR Hypergeom. Function \n');
fprintf('CP-OFDM        |%9.1f dB  |%9.1f dB  |%12.1f dB |%12.1f dB  \n', SIR_OFDM_Matrix       , SIR_OFDM_Object       , SIR_OFDM_Simulation       , SIR_OFDM_TheoryJakes      );
fprintf('OFDM (no CPO)  |%9.1f dB  |%9.1f dB  |%12.1f dB | \n'          , SIR_OFDMnoCP_Matrix   , SIR_OFDMnoCP_Object   , SIR_OFDMnoCP_Simulation                               );
fprintf('FBMC-OQAM      |%9.1f dB  |%9.1f dB  |%12.1f dB | \n'          , SIR_FBMC_Matrix       , SIR_FBMC_Object       , SIR_FBMC_Simulation                                   );
fprintf('==============================================================================================\n');





