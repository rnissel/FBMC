% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script simulates the bit error ratio over SNR in a time-variant 
% multipath channel (doubly-selective channel) for different multicarrier
% modulation techniques:
%       1) OFDM with CP (worst spectral behaviour)
%       2) FBMC (best spectral properties, complex orthogonality is replaced by real orthogonality)
%       3) WOLA (windowed OFDM. The windowing is done at TX and RX)
%       4) FOFDM (filtered OFDM, sinc + Hann, filtering at TX and RX)
%       5) UFMC (filtered OFDM, subband-wise filtering, Dolph-Chebyshev window, cyclic prefix and zero padding are supported, filtering at TX and RX)
% Note that WOLA, FOFDM and UFMC have a high degree of freedom! Depending
% on the window/filter length and CP/ZP length they perform differently.

clear; close all;
addpath('./Theory');

% Sampling rate: should approximately match the power delay profile of the 
% channel and must be larger than "SubcarrierSpacing*NumberOfSubcarriers".
SamplingRate                = 15e3*14*14;
dt                          = 1/SamplingRate;

% Simulation parameters
Simulation_SNR_OFDM_dB            = [-5:5:30];                          % SNR for OFDM in dB. The average transmit power of all methods is the same! However, the SNR might be different due to filtering (in FOFDM and UFMC) or because a different bandwidth is used (different subcarrier spacing or different number of subcarriers).
Simulation_MonteCarloRepetitions  = 500;                                % Number of Monte Carlo repetitions over which we take the average                  
Simulation_PlotSignalPowerTheory  = true;                               % If true, plot also the expected transmit power over time (theory)

% Channel parameters
Channel_Velocity_kmh         = 130;                                     % Velocity in km/h. Note that [mph]*1.6=[kmh] and [m/s]*3.6=[kmh]
Channel_PowerDelayProfile    = 'VehicularA';                            % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
Channel_DopplerModel         = 'Jakes';                                 % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This however only works accuratly if the number of samples and the velocity is sufficiently high.                                       
Channel_CarrierFrequency     = 2.5e9;                                   % Carrier Frequency (Hz)

% General modulation paramters
QAM_ModulationOrder          = 64;                                      % QAM sigal constellation order: 4, 16, 64, 256, 1024,...

% FBMC parameters
FBMC_NumberOfSubcarriers     = 24;                                      % Number of subcarriers
FBMC_NumberOfSymbolsInTime   = 30;                                      % Number FBMC symbols in time
FBMC_SubcarrierSpacing       = 15e3;                                    % Subcarrier spacing (Hz)
FBMC_PrototypeFilter         = 'Hermite-OQAM';                          % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM.
FBMC_OverlappingFactor       = 4;                                       % Overlapping factor, 2,3,4,...

% OFDM parameters
OFDM_NumberOfSubcarriers     = 24;                                      % Number of subcarriers 
OFDM_NumberOfSymbolsInTime   = 14;                                      % Number OFDM symbols in time
OFDM_SubcarrierSpacing       = 15e3;                                    % Subcarrier spacing (Hz)
OFDM_CyclicPrefixLength      = 1/(14*OFDM_SubcarrierSpacing);           % Length of the cyclic prefix (s)

% WOLA parameters
WOLA_NumberOfSubcarriers     = 24;                                      % Number of subcarriers                                   
WOLA_NumberOfSymbolsInTime   = 14;                                      % Number WOLA symbols in time  
WOLA_SubcarrierSpacing       = 15e3;                                    % Subcarrier spacing (Hz)
WOLA_CyclicPrefixLength      = 0;                                       % Length of the cyclic prefix (s) to combat the channel
WOLA_WindowLengthTX          = 1/(14*2*WOLA_SubcarrierSpacing);         % Length of the window overlapping (s) at the transmitter 
WOLA_WindowLengthRX          = 1/(14*2*WOLA_SubcarrierSpacing);         % Length of the window overlapping (s) at the receiver

% FOFDM parameters
FOFDM_NumberOfSubcarriers     = 24;                                     % Number of subcarriers
FOFDM_NumberOfSymbolsInTime   = 14;                                     % Number FOFDM symbols in time                        
FOFDM_SubcarrierSpacing       = 15e3;                                   % Subcarrier spacing (Hz)
FOFDM_CyclicPrefixLength      = 0;                                      % Length of the cyclic prefix (s) to combat the channel
FOFDM_FilterLengthTX          = 0.2*1/(FOFDM_SubcarrierSpacing);        % Length of the transmit filter (s)
FOFDM_FilterLengthRX          = 0.2*1/(FOFDM_SubcarrierSpacing);        % Length of the receive filter (s) 
FOFDM_FilterCylicPrefixLength = 1/(14*FOFDM_SubcarrierSpacing);         % Length of the additional cyclic prefix (s) to combat ISI and ICI due to the filtering

% UFMC parameters
UFMC_NumberOfSubcarriers     = 24;                                      % Number of subcarriers
UFMC_NumberOfSymbolsInTime   = 14;                                      % Number UFMC symbols in time
UFMC_SubcarrierSpacing       = 15e3;                                    % Subcarrier spacing (Hz)
UFMC_CyclicPrefixLength      = 0;                                       % Length of the cyclic prefix (s) to combat the channel. If zero padding is used, this length reprents the zero guard length instead of the CP length.
UFMC_FilterLengthTX          = 1/14*1/(UFMC_SubcarrierSpacing);         % Length of the transmit filter (s)
UFMC_FilterLengthRX          = 1/14*1/(UFMC_SubcarrierSpacing);         % Length of the receive filter (s)
UFMC_FilterCylicPrefixLength = 1/(14*UFMC_SubcarrierSpacing);           % Length of the additional cyclic prefix (or zero guard symbol if ZP is used) in seconds (s). Needed to combat ISI and ICI due to the filtering. However, small ICI and ISI is perfectly feasibly.
UFMC_ZeroPaddingInsteadOfCP  = true;                                    % TRUE for Zero Padding (ZP) and FALSE for a conventional Cyclic Prefix (CP)



%% Generate " +Modulation\" Objects
% FBMC Object
FBMC = Modulation.FBMC(...
    FBMC_NumberOfSubcarriers,...                                        % Number of subcarriers
    FBMC_NumberOfSymbolsInTime,...                                      % Number FBMC symbols in time
    FBMC_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz).  Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    FBMC_PrototypeFilter,...                                            % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM. The data rate of QAM is reduced by a factor of two compared to OQAM, but robustness in doubly-selective channels is inceased
    FBMC_OverlappingFactor, ...                                         % Overlapping factor (also determines oversampling in the frequency domain)                                   
    0, ...                                                              % Initial phase shift
    true ...                                                            % Polyphase implementation
    );
FBMC_BlockOverlapTime = (FBMC.PrototypeFilter.OverlappingFactor-1/2)*FBMC.PHY.TimeSpacing;

% OFDM Object
OFDM = Modulation.OFDM(...
    OFDM_NumberOfSubcarriers,...                                        % Number of subcarriers
    OFDM_NumberOfSymbolsInTime,...                                      % Number OFDM symbols in time                                                 
    OFDM_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz) 
    SamplingRate,...                                                    % Sampling rate (Samples/s)                                       
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    OFDM_CyclicPrefixLength, ...                                        % Length of the cyclic prefix (s)                 
    FBMC_BlockOverlapTime ...                                           % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    );

% Windowed OFDM (WOLA)
WOLA = Modulation.WOLA(...
    WOLA_NumberOfSubcarriers,...                                        % Number subcarriers
    WOLA_NumberOfSymbolsInTime,...                                      % Number WOLA symbols in time 
    WOLA_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    0, ...                                                              % Length of the cyclic prefix (s)
    FBMC_BlockOverlapTime-WOLA_WindowLengthTX/2, ...                    % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    WOLA_WindowLengthTX, ...                                            % Length of the window overlapping (s) at the transmitter 
    WOLA_WindowLengthRX ...                                             % Length of the window overlapping (s) at the receiver
    );

% FOFDM (Filtered OFDM)
FOFDM = Modulation.FOFDM(...
    FOFDM_NumberOfSubcarriers,...                                       % Number of subcarriers
    FOFDM_NumberOfSymbolsInTime,...                                     % Number FOFDM symbols in time                
    FOFDM_SubcarrierSpacing,...                                         % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    0, ...                                                              % Length of the cyclic prefix (s)
    FBMC_BlockOverlapTime-FOFDM_FilterLengthTX/2, ...                   % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission                    
    FOFDM_FilterLengthTX, ...                                           % Length of the transmit filter (s)
    FOFDM_FilterLengthRX, ...                                           % Length of the receive filter (s) 
    FOFDM_FilterCylicPrefixLength ...                                   % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
);

% UFMC (Subband Filtered OFDM)
UFMC = Modulation.UFMC(...
    UFMC_NumberOfSubcarriers,...                                        % Number of subcarriers
    UFMC_NumberOfSymbolsInTime,...                                      % Number UFMC symbols in time
    UFMC_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    0, ...                                                              % Length of the cyclic prefix (s). If zero padding is used, this length reprents the zero guard length instead of the CP length
    FBMC_BlockOverlapTime-UFMC_FilterLengthTX/2, ...                    % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    UFMC_FilterLengthTX, ...                                            % Length of the transmit filter (s)
    UFMC_FilterLengthRX, ...                                            % Length of the receive filter (s)
    UFMC_FilterCylicPrefixLength, ...                                   % Length of the additional cyclic prefix (or zero guard symbol if ZP is used) in seconds (s). Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
    UFMC_ZeroPaddingInsteadOfCP ...                                     % TRUE for Zero Padding (ZP) and FALSE for a conventional Cyclic Prefix (CP)
);

% Number of samples
N_FBMC  = FBMC.Nr.SamplesTotal;
N_OFDM  = OFDM.Nr.SamplesTotal;
N_WOLA  = WOLA.Nr.SamplesTotal;
N_FOFDM = FOFDM.Nr.SamplesTotal;
N_UFMC  = UFMC.Nr.SamplesTotal;
N       = max([N_FBMC N_OFDM N_WOLA N_FOFDM N_UFMC]);

ChannelModel = Channel.FastFading(...
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    Channel_PowerDelayProfile,...                                       % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
    N,...                                                               % Number of total samples
    Channel_Velocity_kmh/3.6*Channel_CarrierFrequency/2.998e8,...       % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    Channel_DopplerModel,...                                            % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
    200, ...                                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
    1,...                                                               % Number of transmit antennas
    1,...                                                               % Number of receive antennas
    true ...                                                            % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );

% PAM and QAM Object
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
if strcmp(FBMC.Method(end-3),'O')
    % FBMC-OQAM transmission, only real valued data symbols
    PAMorQAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
else
    % FBMC-QAM transmission,  complex valued data symbols
    PAMorQAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
end


% Pre-allocate Transmit Power
Ps_FBMC     = zeros(N_FBMC,1);
Ps_OFDM     = zeros(N_OFDM,1);
Ps_WOLA     = zeros(N_WOLA,1);
Ps_FOFDM    = zeros(N_FOFDM,1);
Ps_UFMC     = zeros(N_UFMC,1);

% Pre-allocate Power Spectral Density
PSD_FBMC    = zeros(N_FBMC,1);
PSD_OFDM    = zeros(N_OFDM,1);
PSD_WOLA    = zeros(N_WOLA,1);
PSD_FOFDM   = zeros(N_FOFDM,1);
PSD_UFMC    = zeros(N_UFMC,1);

%% Start Simulation
tic
for i_rep = 1:Simulation_MonteCarloRepetitions
    % Update channel
    ChannelModel.NewRealization;
    
    % Binary data
    BinaryDataStream_FBMC   = randi([0 1], FBMC.Nr.Subcarriers  * FBMC.Nr.MCSymbols  * log2(PAMorQAM.ModulationOrder),1);
    BinaryDataStream_OFDM   = randi([0 1], OFDM.Nr.Subcarriers  * OFDM.Nr.MCSymbols  * log2(QAM.ModulationOrder),1);
    BinaryDataStream_WOLA   = randi([0 1], WOLA.Nr.Subcarriers  * WOLA.Nr.MCSymbols  * log2(QAM.ModulationOrder),1);
    BinaryDataStream_FOFDM  = randi([0 1], FOFDM.Nr.Subcarriers * FOFDM.Nr.MCSymbols * log2(QAM.ModulationOrder),1);
    BinaryDataStream_UFMC   = randi([0 1], FOFDM.Nr.Subcarriers * FOFDM.Nr.MCSymbols * log2(QAM.ModulationOrder),1);
     
    % Transmitted data symbols (map binary data to symbol)
    x_FBMC  = reshape( PAMorQAM.Bit2Symbol(BinaryDataStream_FBMC)  , FBMC.Nr.Subcarriers  , FBMC.Nr.MCSymbols);
    x_OFDM  = reshape(      QAM.Bit2Symbol(BinaryDataStream_OFDM)  , OFDM.Nr.Subcarriers  , OFDM.Nr.MCSymbols);
    x_WOLA  = reshape(      QAM.Bit2Symbol(BinaryDataStream_WOLA)  , WOLA.Nr.Subcarriers  , WOLA.Nr.MCSymbols);
    x_FOFDM = reshape(      QAM.Bit2Symbol(BinaryDataStream_FOFDM) , FOFDM.Nr.Subcarriers , FOFDM.Nr.MCSymbols);
    x_UFMC  = reshape(      QAM.Bit2Symbol(BinaryDataStream_UFMC)  , UFMC.Nr.Subcarriers  , UFMC.Nr.MCSymbols);
        
    % Transmitted signal in the time domain
    s_FBMC  =  FBMC.Modulation( x_FBMC );
    s_OFDM  =  OFDM.Modulation( x_OFDM );
    s_WOLA  =  WOLA.Modulation( x_WOLA );
    s_FOFDM = FOFDM.Modulation( x_FOFDM );
    s_UFMC  =  UFMC.Modulation( x_UFMC );
  
    % Channel convolution
    r_FBMC_noNoise  = ChannelModel.Convolution( s_FBMC );
    r_OFDM_noNoise  = ChannelModel.Convolution( s_OFDM );
    r_WOLA_noNoise  = ChannelModel.Convolution( s_WOLA );
    r_FOFDM_noNoise = ChannelModel.Convolution( s_FOFDM );
    r_UFMC_noNoise  = ChannelModel.Convolution( s_UFMC );

    % Channel for one tap equalizer (represents only an approximation because pulseform is not taken into account)
    h_FBMC  = ChannelModel.GetTransferFunction( FBMC.GetTimeIndexMidPos  , FBMC.Implementation.FFTSize  , (1:FBMC.Nr.Subcarriers)  + FBMC.Implementation.IntermediateFrequency);
    h_OFDM  = ChannelModel.GetTransferFunction( OFDM.GetTimeIndexMidPos  , OFDM.Implementation.FFTSize  , (1:OFDM.Nr.Subcarriers)  + OFDM.Implementation.IntermediateFrequency);
    h_WOLA  = ChannelModel.GetTransferFunction( WOLA.GetTimeIndexMidPos  , WOLA.Implementation.FFTSize  , (1:WOLA.Nr.Subcarriers)  + WOLA.Implementation.IntermediateFrequency);
    h_FOFDM = ChannelModel.GetTransferFunction( FOFDM.GetTimeIndexMidPos , FOFDM.Implementation.FFTSize , (1:FOFDM.Nr.Subcarriers) + FOFDM.Implementation.IntermediateFrequency);
    h_UFMC  = ChannelModel.GetTransferFunction( UFMC.GetTimeIndexMidPos  , UFMC.Implementation.FFTSize  , (1:UFMC.Nr.Subcarriers)  + UFMC.Implementation.IntermediateFrequency);

    % Calculate the transmitted power over time
    Ps_FBMC  = Ps_FBMC  + abs(s_FBMC).^2;
    Ps_OFDM  = Ps_OFDM  + abs(s_OFDM).^2;
    Ps_WOLA  = Ps_WOLA  + abs(s_WOLA).^2;
    Ps_FOFDM = Ps_FOFDM + abs(s_FOFDM).^2;
    Ps_UFMC  = Ps_UFMC  + abs(s_UFMC).^2;
    
    % Calculat the power spectral density
    PSD_FBMC  = PSD_FBMC  + abs(fft(s_FBMC)/sqrt(N_FBMC)).^2;
    PSD_OFDM  = PSD_OFDM  + abs(fft(s_OFDM)/sqrt(N_OFDM)).^2;
    PSD_WOLA  = PSD_WOLA  + abs(fft(s_WOLA)/sqrt(N_WOLA)).^2;
    PSD_FOFDM = PSD_FOFDM + abs(fft(s_FOFDM)/sqrt(N_FOFDM)).^2;
    PSD_UFMC  = PSD_UFMC  + abs(fft(s_UFMC)/sqrt(N_UFMC)).^2;    
    
for i_SNR = 1:length(Simulation_SNR_OFDM_dB)
    % Add noise
    SNR_OFDM_dB = Simulation_SNR_OFDM_dB(i_SNR);
    Pn_time     = 1/OFDM.GetSymbolNoisePower(1)*10^(-SNR_OFDM_dB/10);
    noise       = sqrt(Pn_time/2)*(randn(N,1)+1j*randn(N,1));
        
    r_FBMC  = r_FBMC_noNoise  + noise(1:N_FBMC);
    r_OFDM  = r_OFDM_noNoise  + noise(1:N_OFDM);
    r_WOLA  = r_WOLA_noNoise  + noise(1:N_WOLA);
    r_FOFDM = r_FOFDM_noNoise + noise(1:N_FOFDM);
    r_UFMC  = r_UFMC_noNoise  + noise(1:N_UFMC);
   
    % Demodulate FBMC signal
    y_FBMC  =  FBMC.Demodulation(r_FBMC);
    y_OFDM  =  OFDM.Demodulation(r_OFDM);
    y_WOLA  =  WOLA.Demodulation(r_WOLA);
    y_FOFDM = FOFDM.Demodulation(r_FOFDM);
    y_UFMC  =  UFMC.Demodulation(r_UFMC);
    
    % One-tap equalizer
    y_OneTapEqualizer_FBMC  =  y_FBMC./h_FBMC;
    y_OneTapEqualizer_OFDM  =  y_OFDM./h_OFDM;
    y_OneTapEqualizer_WOLA  =  y_WOLA./h_WOLA;
    y_OneTapEqualizer_FOFDM = y_FOFDM./h_FOFDM;
    y_OneTapEqualizer_UFMC  =  y_UFMC./h_UFMC;
    
    % Detect symbols (quantization and demapping to bits)
    DetectedBitStream_OneTapEqualizer_FBMC  = PAMorQAM.Symbol2Bit(y_OneTapEqualizer_FBMC(:));
    DetectedBitStream_OneTapEqualizer_OFDM  = QAM.Symbol2Bit(y_OneTapEqualizer_OFDM(:));
    DetectedBitStream_OneTapEqualizer_WOLA  = QAM.Symbol2Bit(y_OneTapEqualizer_WOLA(:));
    DetectedBitStream_OneTapEqualizer_FOFDM = QAM.Symbol2Bit(y_OneTapEqualizer_FOFDM(:));
    DetectedBitStream_OneTapEqualizer_UFMC  = QAM.Symbol2Bit(y_OneTapEqualizer_UFMC(:));
   
    % Calculate the BER
    BER_FBMC_OneTapEqualizer(i_SNR,i_rep)   = mean( BinaryDataStream_FBMC~=DetectedBitStream_OneTapEqualizer_FBMC );
    BER_OFDM_OneTapEqualizer(i_SNR,i_rep)   = mean( BinaryDataStream_OFDM~=DetectedBitStream_OneTapEqualizer_OFDM );
    BER_WOLA_OneTapEqualizer(i_SNR,i_rep)   = mean( BinaryDataStream_WOLA~=DetectedBitStream_OneTapEqualizer_WOLA );
    BER_FOFDM_OneTapEqualizer(i_SNR,i_rep)  = mean( BinaryDataStream_FOFDM~=DetectedBitStream_OneTapEqualizer_FOFDM );
    BER_UFMC_OneTapEqualizer(i_SNR,i_rep)   = mean( BinaryDataStream_UFMC~=DetectedBitStream_OneTapEqualizer_UFMC );  

end
TimeNeededSoFar = toc;
if mod(i_rep,25)==0
    disp([int2str(i_rep/Simulation_MonteCarloRepetitions*100) '% Completed. Time Left, approx. ' int2str(TimeNeededSoFar/i_rep*(Simulation_MonteCarloRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(Simulation_MonteCarloRepetitions-i_rep)/3600) 'hour'])
end
end

% Take "average"
Ps_FBMC  = Ps_FBMC/Simulation_MonteCarloRepetitions;
Ps_OFDM  = Ps_OFDM/Simulation_MonteCarloRepetitions;
Ps_WOLA  = Ps_WOLA/Simulation_MonteCarloRepetitions;
Ps_FOFDM = Ps_FOFDM/Simulation_MonteCarloRepetitions;
Ps_UFMC  = Ps_UFMC/Simulation_MonteCarloRepetitions;
PSD_FBMC  = PSD_FBMC/Simulation_MonteCarloRepetitions;
PSD_OFDM  = PSD_OFDM/Simulation_MonteCarloRepetitions;
PSD_WOLA  = PSD_WOLA/Simulation_MonteCarloRepetitions;
PSD_FOFDM = PSD_FOFDM/Simulation_MonteCarloRepetitions;
PSD_UFMC  = PSD_UFMC/Simulation_MonteCarloRepetitions; 
    

% Define colors for different modulation schemes
ColorFBMC   = [0 0 1]*0.5;
ColorOFDM   = [1 0 0];
ColorWOLA   = [1 0 1];
ColorFOFDM  = [1 1 0]*0.7;
ColorUFMC   = [0 1 0]*0.5;

% Plot BER
figure(1);
semilogy( Simulation_SNR_OFDM_dB , mean(BER_FBMC_OneTapEqualizer,2)  ,'x-','color',ColorFBMC);  hold on;
semilogy( Simulation_SNR_OFDM_dB , mean(BER_OFDM_OneTapEqualizer,2)  ,'o-','color',ColorOFDM);  hold on;
semilogy( Simulation_SNR_OFDM_dB , mean(BER_WOLA_OneTapEqualizer,2)  ,'s-','color',ColorWOLA);  hold on;
semilogy( Simulation_SNR_OFDM_dB , mean(BER_FOFDM_OneTapEqualizer,2) ,'d-','color',ColorFOFDM); hold on;
semilogy( Simulation_SNR_OFDM_dB , mean(BER_UFMC_OneTapEqualizer,2)  ,'*-','color',ColorUFMC);  hold on;
ylabel('Bit Error Ratio');
xlabel('Signal-to-Noise Ratio for OFDM [dB]');
% Theoretical bit error probability
if strcmp(Channel_PowerDelayProfile,'AWGN')
    SNR_morePoints = linspace(min(Simulation_SNR_OFDM_dB),max(Simulation_SNR_OFDM_dB),100);
    BitErrorProability = BitErrorProbabilityAWGN(SNR_morePoints,QAM.SymbolMapping,QAM.BitMapping);
    semilogy( SNR_morePoints , BitErrorProability  ,'black');  hold on;
    ylim([10^-5 1]);
    legend({'FBMC','OFDM','WOLA','f-OFDM','UFMC','Theory(OFDM)'},'Location','SouthWest');
else
    SNR_morePoints = linspace(min(Simulation_SNR_OFDM_dB),max(Simulation_SNR_OFDM_dB),100);    
    BitErrorProability = BitErrorProbabilityDoublyFlatRayleigh(SNR_morePoints,QAM.SymbolMapping,QAM.BitMapping);
    semilogy( SNR_morePoints , BitErrorProability  ,'black');  hold on;
    legend({'FBMC','OFDM','WOLA','f-OFDM','UFMC','Theory (Doubly-Flat, OFDM)'},'Location','SouthWest');    
end

% Plot simulated power
figure(2);
time_FBMC   = ((0:N_FBMC-1) *dt-FBMC_BlockOverlapTime);
time_OFDM   = ((0:N_OFDM-1) *dt-FBMC_BlockOverlapTime);
time_WOLA   = ((0:N_WOLA-1) *dt-FBMC_BlockOverlapTime);
time_FOFDM  = ((0:N_FOFDM-1)*dt-FBMC_BlockOverlapTime);
time_UFMC   = ((0:N_UFMC-1) *dt-FBMC_BlockOverlapTime);
plot( time_FBMC/1e-3  , Ps_FBMC,  'color',ColorFBMC); hold on;
plot( time_OFDM/1e-3  , Ps_OFDM,  'color',ColorOFDM); hold on;
plot( time_WOLA/1e-3  , Ps_WOLA,  'color',ColorWOLA); hold on;
plot( time_FOFDM/1e-3 , Ps_FOFDM, 'color',ColorFOFDM); hold on;
plot( time_UFMC/1e-3  , Ps_UFMC,  'color',ColorUFMC); hold on;
ylabel('Signal Power');
xlabel('Time [ms]');
legend({'FBMC','OFDM','WOLA','f-OFDM','UFMC'});
title('Simulation');
% Theory
if Simulation_PlotSignalPowerTheory
    figure(5);
    plot( time_FBMC/1e-3  , FBMC.PlotTransmitPower ,  'color',ColorFBMC); hold on;
    plot( time_OFDM/1e-3  , OFDM.PlotTransmitPower ,  'color',ColorOFDM); hold on;
    plot( time_WOLA/1e-3  , WOLA.PlotTransmitPower ,  'color',ColorWOLA); hold on;
    plot( time_FOFDM/1e-3 , FOFDM.PlotTransmitPower,  'color',ColorFOFDM); hold on;
    plot( time_UFMC/1e-3  , UFMC.PlotTransmitPower ,  'color',ColorUFMC); hold on;
    ylabel('Signal Power');
    xlabel('Time [ms]');
    legend({'FBMC','OFDM','WOLA','f-OFDM','UFMC'});
    title('Theory');
end


figure(3);
frequency_FBMC  = (0:N_FBMC  -1) / (N_FBMC *dt);
frequency_OFDM  = (0:N_OFDM  -1) / (N_OFDM *dt);
frequency_WOLA  = (0:N_WOLA  -1) / (N_WOLA *dt);
frequency_FOFDM = (0:N_FOFDM -1) / (N_FOFDM*dt);
frequency_UFMC  = (0:N_UFMC  -1) / (N_UFMC *dt);
plot( frequency_FBMC /1e6 , 10*log10(PSD_FBMC/max(PSD_FBMC)) ,'color',ColorFBMC);  hold on;
plot( frequency_OFDM /1e6 , 10*log10(PSD_OFDM/max(PSD_FBMC)) ,'color',ColorOFDM);  hold on;
plot( frequency_WOLA /1e6 , 10*log10(PSD_WOLA/max(PSD_FBMC)) ,'color',ColorWOLA);  hold on;
plot( frequency_FOFDM/1e6 , 10*log10(PSD_FOFDM/max(PSD_FBMC)),'color',ColorFOFDM); hold on;
plot( frequency_UFMC /1e6 , 10*log10(PSD_UFMC/max(PSD_FBMC)) ,'color',ColorUFMC);  hold on;
ylabel('Power Spectral Density [dB]');
xlabel('Frequency [MHz]');
ylim([-100 0]);
legend({'FBMC','OFDM','WOLA','f-OFDM','UFMC'});
title('Simulation');


% Note that the bitrate assumes infinitly many time symbols
fprintf('===========================================================================================\n');
fprintf('================== Basic Settings (guard time and band are ignored) =======================\n');
fprintf('               |(complex)TF-Spacing| Bandwidth(FL)|    Time(KT) |      BitRate \n');
fprintf('OFDM (with CP) |%17.2f  |%8.2f MHz  |%8.2f ms  |%8.2f Mbit/s  | \n', OFDM.PHY.TimeSpacing*OFDM.PHY.SubcarrierSpacing   , OFDM.PHY.SubcarrierSpacing*OFDM.Nr.Subcarriers/1e6   , OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols/1e-3   , length(DetectedBitStream_OneTapEqualizer_OFDM) /(OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols)/1e6   );
fprintf('FBMC-OQAM      |%17.2f  |%8.2f MHz  |%8.2f ms  |%8.2f Mbit/s  | \n', FBMC.PHY.TimeSpacing*FBMC.PHY.SubcarrierSpacing*2 , FBMC.PHY.SubcarrierSpacing*FBMC.Nr.Subcarriers/1e6   , FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols/1e-3   , length(DetectedBitStream_OneTapEqualizer_FBMC) /(FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols)/1e6   );
fprintf('WOLA           |%17.2f  |%8.2f MHz  |%8.2f ms  |%8.2f Mbit/s  | \n', WOLA.PHY.TimeSpacing*WOLA.PHY.SubcarrierSpacing   , WOLA.PHY.SubcarrierSpacing*WOLA.Nr.Subcarriers/1e6   , WOLA.PHY.TimeSpacing*WOLA.Nr.MCSymbols/1e-3   , length(DetectedBitStream_OneTapEqualizer_WOLA) /(WOLA.PHY.TimeSpacing*WOLA.Nr.MCSymbols)/1e6   );
fprintf('FOFDM          |%17.2f  |%8.2f MHz  |%8.2f ms  |%8.2f Mbit/s  | \n', FOFDM.PHY.TimeSpacing*FOFDM.PHY.SubcarrierSpacing , FOFDM.PHY.SubcarrierSpacing*FOFDM.Nr.Subcarriers/1e6 , FOFDM.PHY.TimeSpacing*FOFDM.Nr.MCSymbols/1e-3 , length(DetectedBitStream_OneTapEqualizer_FOFDM)/(FOFDM.PHY.TimeSpacing*FOFDM.Nr.MCSymbols)/1e6 );
fprintf('UFMC           |%17.2f  |%8.2f MHz  |%8.2f ms  |%8.2f Mbit/s  | \n', UFMC.PHY.TimeSpacing*UFMC.PHY.SubcarrierSpacing   , UFMC.PHY.SubcarrierSpacing*UFMC.Nr.Subcarriers/1e6   , UFMC.PHY.TimeSpacing*UFMC.Nr.MCSymbols/1e-3   , length(DetectedBitStream_OneTapEqualizer_UFMC) /(UFMC.PHY.TimeSpacing*UFMC.Nr.MCSymbols)/1e6   );
fprintf('===========================================================================================\n');
fprintf('===========================================================================================\n');


fprintf('=============================== Transmitted Power =========================================\n');
fprintf('               | (Sim.) Trans. Energy |     Time    | (Sim.) Av. Power | SNR rel. to OFDM | \n');
fprintf('OFDM (with CP) | %17.2f mJ |%8.2f ms  |%14.2f W  |%14.2f dB | \n', sum(Ps_OFDM)*dt/1e-3   , OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols/1e-3   , sum(Ps_OFDM)*dt/(OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols)    , 10*log10(OFDM.GetSymbolNoisePower(1) / OFDM.GetSymbolNoisePower(1))  );
fprintf('FBMC-OQAM      | %17.2f mJ |%8.2f ms  |%14.2f W  |%14.2f dB | \n', sum(Ps_FBMC)*dt/1e-3   , FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols/1e-3   , sum(Ps_FBMC)*dt/(FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols)    , 10*log10(OFDM.GetSymbolNoisePower(1) / FBMC.GetSymbolNoisePower(1)*(1+strcmp(FBMC.Method(end-3),'O'))));
fprintf('WOLA           | %17.2f mJ |%8.2f ms  |%14.2f W  |%14.2f dB | \n', sum(Ps_WOLA)*dt/1e-3   , WOLA.PHY.TimeSpacing*WOLA.Nr.MCSymbols/1e-3   , sum(Ps_WOLA)*dt/(WOLA.PHY.TimeSpacing*WOLA.Nr.MCSymbols)    , 10*log10(OFDM.GetSymbolNoisePower(1) / WOLA.GetSymbolNoisePower(1))  );
fprintf('FOFDM          | %17.2f mJ |%8.2f ms  |%14.2f W  |%14.2f dB | \n', sum(Ps_FOFDM)*dt/1e-3  , FOFDM.PHY.TimeSpacing*FOFDM.Nr.MCSymbols/1e-3 , sum(Ps_FOFDM)*dt/(FOFDM.PHY.TimeSpacing*FOFDM.Nr.MCSymbols) , 10*log10(OFDM.GetSymbolNoisePower(1) / mean(mean(FOFDM.GetSymbolNoisePower(1))) )  );
fprintf('UFMC           | %17.2f mJ |%8.2f ms  |%14.2f W  |%14.2f dB | \n', sum(Ps_UFMC)*dt/1e-3   , UFMC.PHY.TimeSpacing*UFMC.Nr.MCSymbols/1e-3   , sum(Ps_UFMC)*dt/(UFMC.PHY.TimeSpacing*UFMC.Nr.MCSymbols)    , 10*log10(OFDM.GetSymbolNoisePower(1) / mean(mean(UFMC.GetSymbolNoisePower(1))) )  );
fprintf('===========================================================================================\n');
fprintf('===========================================================================================\n');


figure(4);
% Plot also the power delay profile.
ChannelModel.PlotPowerDelayProfile;








