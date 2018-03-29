% =========================================================================   
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =========================================================================      
% This script emulates our MIMO BER measurements. MIMO in FBMC-OQAM is 
% enabled by spreading data symbols in time (or frequency). We consider
% three transmit blocks in time but evaluate the performance only for the
% second block => emulates infinitly many blocks. Note that the channel
% estimation is optimized for a low delay spread and a low doppler spread 
% and might not work for different channel parameters. If you want to
% simulate other channel conditions, best to set  
% Simulation_IncludePerfectCSI to true in order to avoid any errors due to
% the channel estimation error.
% More information about MIMO in FBMC-OQAM can be found in 
% "Enabling Low-Complexity MIMO in FBMC-OQAM", R. Nissel, et.al

clear; close all;

%% Parameters
% Simulation parameters
Simulation_MonteCarloRepetitions            = 1024;                     % Number of Monte Carlo repetitions
Simulation_SNR_dB_FBMC                      = -5:2.5:20;                % Signal-to-noise ratio for FBMC
Simulation_IncludePerfectCSI                = false;                    % If set to true, also calculate the BER for perfect channel knowledge (takes more time because we calculate the true one tap channel, using matrices, and not just an approximation!)

% Channel parameters are similar to our measurement at 2.5GHz (Rayleigh 
% fading). For our 60GHz measurements, we observe Rician fading (not
% included here)
Channel_PowerDelayProfile                   = 'PedestrianA';            % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)     
Channel_Velocity_kmh                        = 0;                        % Velocity 
Channel_CarrierFrequency                    = 2.5e9;                    % Carrier frequency (has no influence if the velocity is 0!) 

% Modulation parameters:
L                                           = 12;                       % Number of subcarriers
K                                           = 2^5;                      % Number of FBMC symbols per block = spreading length. Must be a power of two: 2^2, 2^3, 2^4, 2^5,... 
QAM_ModulationOrder                         = 16;                       % Modulation order, 4, 16, 64, 256, 1024,...
SubcarrierSpacing                           = 15e3;                     % Subcarrier spacing (15kHz, same as LTE)


%% Objects
Kall            = K*3+3;                                                % Total number of FBMC symbols (3 guard symbols!)
SamplingRate    =  SubcarrierSpacing*L*K/4;                             % Sampling rate (must be larger than the subcarrier spacing times L)
% FBMC Object
FBMC = Modulation.FBMC(...
    L,...                                                               % Number of subcarriers
    Kall,...                                                            % Number of FBMC symbols
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency first subcarrier (Hz)
    false,...                                                           % Transmit real valued signal
    'Hermite-OQAM',...                                                  % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    8, ...                                                              % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                              % Initial phase shift
    true ...                                                            % Polyphase implementation
    );
% OFDM Object
OFDM = Modulation.OFDM(...
    L,...                                                               % Number of subcarriers
    K/2*3,...                                                           % Number of OFDM Symbols
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency first subcarrier (Hz)
    false,...                                                           % Transmit real valued signal
    1/(SubcarrierSpacing*K), ...                                        % Cyclic prefix length (s)
    (8-1/2)*1/SubcarrierSpacing*1/2 ...                                 % Zero guard length (s)
    );
if FBMC.Nr.SamplesTotal == OFDM.Nr.SamplesTotal
    N = FBMC.Nr.SamplesTotal;
else
    error('Number of samples in OFDM and FBMC have to be the same.');
end
ChannelModel = Channel.FastFading(...
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    Channel_PowerDelayProfile,...                                       % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
    N,...                                                               % Number of total samples
    Channel_Velocity_kmh/3.6*Channel_CarrierFrequency/2.998e8,...       % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    'Jakes',...                                                         % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
    200, ...                                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
    2,...                                                               % Number of transmit antennas
    2,...                                                               % Number of receive antennas
    true ...                                                            % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );

% Alamouti Object
Alamouti = MIMO.SpaceCoding(...
    'Alamouti2x1',...                                                   % Only Alamotui2x1 is currently implemented!
    1 ...                                                               % Frequency spreading = 0; time spreading = 1
    );

% Modulation Object
QAM             = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
% For ML Detection
ML_MapIndex1    = reshape(repmat((1:QAM.ModulationOrder),QAM.ModulationOrder,1),1,QAM.ModulationOrder^2);
ML_MapIndex2    = reshape(repmat((1:QAM.ModulationOrder).',1,QAM.ModulationOrder),1,QAM.ModulationOrder^2);
ML_Mapping      = QAM.SymbolMapping([ML_MapIndex1;ML_MapIndex2]);

%% Get coding matrix for QAM transmission in FBMC-OQAM
FBMC.SetNrMCSymbols(K);
C               = FBMC.GetPrecodingMatrixForQAMinOQAM(1,1);
CBlock          = [zeros(L,L*K/2)  ;   C  ]; 
FBMC.SetNrMCSymbols(Kall);  
% Power Normalization in FBMC so that average transmit power is the same in OFDM and FBMC
FBMCPowerNormalization = sqrt(2*(K+1)/K); % 2: Spreading(FBMC-OQAM structure); (K+1)/K due to zero guard symbol

%% Pilot Matrix: 0=DataSymbol, 1=PilotSymbol, -1=ZeroSymbol;
PilotMatrixBlockAntenna1                    = zeros(L,K/2);
PilotMatrixBlockAntenna1(2:6:end,1:8:end)   = 1; 
PilotMatrixBlockAntenna1(5:6:end,5:8:end)   = 1;
PilotMatrixBlockAntenna1(2:6:end,2:8:end)   = -1; 
PilotMatrixBlockAntenna1(5:6:end,6:8:end)   = -1;
PilotMatrixBlockAntenna2                    = PilotMatrixBlockAntenna1*(-1);
PilotMatrixBlockAntenna(:,:,1)              = PilotMatrixBlockAntenna1;
PilotMatrixBlockAntenna(:,:,2)              = PilotMatrixBlockAntenna2;

NrPilots                                    = sum(PilotMatrixBlockAntenna1(:)==1);
NrDataSymbols                               = sum(PilotMatrixBlockAntenna1(:)==0);
NrTransmittedSymbols                        = length(PilotMatrixBlockAntenna1(:));

%% TX and RX matrices of the second block
IndexBlock2_OFDM    = (1+(K/2)*L) : (L+(2*K/2-1)*L);
G_TX_OFDM           = OFDM.GetTXMatrix; 
G_TX_OFDM           = G_TX_OFDM( : , IndexBlock2_OFDM );
G_RX_OFDM           = OFDM.GetRXMatrix; 
G_RX_OFDM           = G_RX_OFDM( IndexBlock2_OFDM , : );

IndexBlock2_FBMC    = (1+(K+1)*L) : (L+(2*K+1)*L);
G_TX_FBMC_Coding    = FBMC.GetTXMatrix; 
G_TX_FBMC_Coding    = G_TX_FBMC_Coding( : , IndexBlock2_FBMC ) * CBlock;
G_RX_FBMC_Coding    = FBMC.GetRXMatrix; 
G_RX_FBMC_Coding    = CBlock' * G_RX_FBMC_Coding( IndexBlock2_FBMC , : );


tic;
%% Simulate Over Different Channel Realizations
for i_rep = 1:Simulation_MonteCarloRepetitions

    %% Generate Data and Pilots
    % Pilot Symbols: The pilot symbol power is increased by a factor of two (pilots for the other antenna are zero)
    x_PilotAntenna1             = QAM.SymbolMapping(randi(QAM.ModulationOrder,NrPilots,1));
    x_PilotAntenna2             = QAM.SymbolMapping(randi(QAM.ModulationOrder,NrPilots,1));
    x_PilotAntenna1             = x_PilotAntenna1./abs(x_PilotAntenna1)*sqrt(2);
    x_PilotAntenna2             = x_PilotAntenna2./abs(x_PilotAntenna2)*sqrt(2);

    % Binary Data Stream
    BinaryDataStream_Alamouti   = randi([0 1],NrDataSymbols*log2(QAM.ModulationOrder),1);
    BinaryDataStream_SMAntenna1 = randi([0 1],NrDataSymbols*log2(QAM.ModulationOrder),1); %Spatial Multiplexing
    BinaryDataStream_SMAntenna2 = randi([0 1],NrDataSymbols*log2(QAM.ModulationOrder),1); %Spatial Multiplexing

    % Transmitted Alamouti Symbols
    x_Alamouti = nan(L,K/2);
    x_Alamouti(PilotMatrixBlockAntenna1==0)         = QAM.Bit2Symbol(BinaryDataStream_Alamouti);
    x_Alamouti_Coded                                = Alamouti.Encoder(x_Alamouti);
    x_Alamouti_Coded(PilotMatrixBlockAntenna==1)    = [x_PilotAntenna1;x_PilotAntenna2];
    x_Alamouti_Coded(PilotMatrixBlockAntenna==-1)   = 0;
    x_Alamouti_Coded_Antenna1                       = x_Alamouti_Coded(:,:,1);
    x_Alamouti_Coded_Antenna2                       = x_Alamouti_Coded(:,:,2);

    % Transmitted Spatial Multiplexed Symbols
    x_SM_Antenna1                                   = nan(L,K/2);
    x_SM_Antenna1(PilotMatrixBlockAntenna1==0)      = QAM.Bit2Symbol(BinaryDataStream_SMAntenna1);
    x_SM_Antenna1(PilotMatrixBlockAntenna1==1)      = x_PilotAntenna1;
    x_SM_Antenna1(PilotMatrixBlockAntenna1==-1)     = 0;
    x_SM_Antenna2                                   = nan(L,K/2);
    x_SM_Antenna2(PilotMatrixBlockAntenna2==0)      = QAM.Bit2Symbol(BinaryDataStream_SMAntenna2);
    x_SM_Antenna2(PilotMatrixBlockAntenna2==1)      = x_PilotAntenna2;
    x_SM_Antenna2(PilotMatrixBlockAntenna2==-1)     = 0;

    % Data symbols of the the first and second block (chosen randomly to keep it simple)
    x_Block1_Antenna1   = QAM.SymbolMapping(randi(QAM.ModulationOrder,L,K/2,1));
    x_Block1_Antenna2   = QAM.SymbolMapping(randi(QAM.ModulationOrder,L,K/2,1));
    x_Block2_Antenna1   = QAM.SymbolMapping(randi(QAM.ModulationOrder,L,K/2,1));
    x_Block2_Antenna2   = QAM.SymbolMapping(randi(QAM.ModulationOrder,L,K/2,1));

    % Account for pilots in block 1 and 3
    x_Block1_Antenna1( PilotMatrixBlockAntenna1== 1 )   = x_Block1_Antenna1(PilotMatrixBlockAntenna1==1)./abs(x_Block1_Antenna1(PilotMatrixBlockAntenna1==1))*sqrt(2);
    x_Block1_Antenna2( PilotMatrixBlockAntenna2== 1 )   = x_Block1_Antenna2(PilotMatrixBlockAntenna2==1)./abs(x_Block1_Antenna2(PilotMatrixBlockAntenna2==1))*sqrt(2);
    x_Block1_Antenna1( PilotMatrixBlockAntenna1==-1 )   = 0;
    x_Block1_Antenna2( PilotMatrixBlockAntenna2==-1 )   = 0;
    x_Block2_Antenna1( PilotMatrixBlockAntenna1== 1 )   = x_Block2_Antenna1(PilotMatrixBlockAntenna1==1)./abs(x_Block2_Antenna1(PilotMatrixBlockAntenna1==1))*sqrt(2);
    x_Block2_Antenna2( PilotMatrixBlockAntenna2== 1 )   = x_Block2_Antenna2(PilotMatrixBlockAntenna2==1)./abs(x_Block2_Antenna2(PilotMatrixBlockAntenna2==1))*sqrt(2);
    x_Block2_Antenna1( PilotMatrixBlockAntenna1==-1 )   = 0;
    x_Block2_Antenna2( PilotMatrixBlockAntenna2==-1 )   = 0;


    %% Transmitted Signal
    TransmittedSymbols_OFDM_Alamouti_Antenna1   = [ x_Block1_Antenna1 x_Alamouti_Coded_Antenna1 x_Block2_Antenna1 ];
    TransmittedSymbols_OFDM_Alamouti_Antenna2   = [ x_Block1_Antenna2 x_Alamouti_Coded_Antenna2 x_Block2_Antenna2 ];
    s_OFDM_Alamouti_Antenna1                    = 1/sqrt(2) * OFDM.Modulation( TransmittedSymbols_OFDM_Alamouti_Antenna1 );
    s_OFDM_Alamouti_Antenna2                    = 1/sqrt(2) * OFDM.Modulation( TransmittedSymbols_OFDM_Alamouti_Antenna2 );

    TransmittedSymbols_FBMC_Alamouti_Antenna1   = reshape(CBlock * [ x_Block1_Antenna1(:) x_Alamouti_Coded_Antenna1(:) x_Block2_Antenna1(:)] , L , K*3+3 );
    TransmittedSymbols_FBMC_Alamouti_Antenna2   = reshape(CBlock * [ x_Block1_Antenna2(:) x_Alamouti_Coded_Antenna2(:) x_Block2_Antenna2(:)] , L , K*3+3 );
    s_FBMC_Alamouti_Antenna1                    = FBMCPowerNormalization/sqrt(2) * FBMC.Modulation( TransmittedSymbols_FBMC_Alamouti_Antenna1 );
    s_FBMC_Alamouti_Antenna2                    = FBMCPowerNormalization/sqrt(2) * FBMC.Modulation( TransmittedSymbols_FBMC_Alamouti_Antenna2 );

    TransmittedSymbols_OFDM_SM_Antenna1         = [ x_Block1_Antenna1 x_SM_Antenna1 x_Block2_Antenna1 ];
    TransmittedSymbols_OFDM_SM_Antenna2         = [ x_Block1_Antenna2 x_SM_Antenna2 x_Block2_Antenna2 ];
    s_OFDM_SM_Antenna1                          = 1/sqrt(2) * OFDM.Modulation( TransmittedSymbols_OFDM_SM_Antenna1 );
    s_OFDM_SM_Antenna2                          = 1/sqrt(2) * OFDM.Modulation( TransmittedSymbols_OFDM_SM_Antenna2 );

    TransmittedSymbols_FBMC_SM_Antenna1         = reshape(CBlock * [ x_Block1_Antenna1(:) x_SM_Antenna1(:) x_Block2_Antenna1(:) ], L , K*3+3 );
    TransmittedSymbols_FBMC_SM_Antenna2         = reshape(CBlock * [ x_Block1_Antenna2(:) x_SM_Antenna2(:) x_Block2_Antenna2(:) ], L , K*3+3 );
    s_FBMC_SM_Antenna1                          = FBMCPowerNormalization/sqrt(2) * FBMC.Modulation( TransmittedSymbols_FBMC_SM_Antenna1 );
    s_FBMC_SM_Antenna2                          = FBMCPowerNormalization/sqrt(2) * FBMC.Modulation( TransmittedSymbols_FBMC_SM_Antenna2 );


    %% Channel
    % Update channel: new realization
    ChannelModel.NewRealization;

    % Received signal without noise
    r_OFDM_Alamouti_noNoise     = ChannelModel.Convolution([s_OFDM_Alamouti_Antenna1 s_OFDM_Alamouti_Antenna2]);
    r_FBMC_Alamouti_noNoise     = ChannelModel.Convolution([s_FBMC_Alamouti_Antenna1 s_FBMC_Alamouti_Antenna2]);
    r_OFDM_SM_noNoise           = ChannelModel.Convolution([s_OFDM_SM_Antenna1 s_OFDM_SM_Antenna2]);
    r_FBMC_SM_noNoise           = ChannelModel.Convolution([s_FBMC_SM_Antenna1 s_FBMC_SM_Antenna2]);

    if Simulation_IncludePerfectCSI
        % Perfect CSI (one-tap) OFDM
        H11_OFDM = reshape(sum((G_RX_OFDM * ChannelModel.GetConvolutionMatrix{1,1}) .* G_TX_OFDM.' ,2 ) ,L , K/2 )/sqrt(2);
        H12_OFDM = reshape(sum((G_RX_OFDM * ChannelModel.GetConvolutionMatrix{1,2}) .* G_TX_OFDM.' ,2 ) ,L , K/2 )/sqrt(2);
        H21_OFDM = reshape(sum((G_RX_OFDM * ChannelModel.GetConvolutionMatrix{2,1}) .* G_TX_OFDM.' ,2 ) ,L , K/2 )/sqrt(2);
        H22_OFDM = reshape(sum((G_RX_OFDM * ChannelModel.GetConvolutionMatrix{2,2}) .* G_TX_OFDM.' ,2 ) ,L , K/2 )/sqrt(2);

        PerfectChannel_OFDM_Alamouti(:,:,1) = H11_OFDM;
        PerfectChannel_OFDM_Alamouti(:,:,2) = H12_OFDM;


        % Perfect CSI (one-tap) FBMC
        H11_FBMC = reshape(sum((G_RX_FBMC_Coding * ChannelModel.GetConvolutionMatrix{1,1}) .* G_TX_FBMC_Coding.' ,2 ) ,L , K/2 ) /sqrt(2) * FBMCPowerNormalization;
        H12_FBMC = reshape(sum((G_RX_FBMC_Coding * ChannelModel.GetConvolutionMatrix{1,2}) .* G_TX_FBMC_Coding.' ,2 ) ,L , K/2 ) /sqrt(2) * FBMCPowerNormalization;
        H21_FBMC = reshape(sum((G_RX_FBMC_Coding * ChannelModel.GetConvolutionMatrix{2,1}) .* G_TX_FBMC_Coding.' ,2 ) ,L , K/2 ) /sqrt(2) * FBMCPowerNormalization;
        H22_FBMC = reshape(sum((G_RX_FBMC_Coding * ChannelModel.GetConvolutionMatrix{2,2}) .* G_TX_FBMC_Coding.' ,2 ) ,L , K/2 ) /sqrt(2) * FBMCPowerNormalization;

        PerfectChannel_FBMC_Alamouti(:,:,1) = H11_FBMC;
        PerfectChannel_FBMC_Alamouti(:,:,2) = H12_FBMC;
    end

    %% Simulate Over Different Noise Values
    for i_SNR = 1:length(Simulation_SNR_dB_FBMC)
        SNR_dB_FBMC                     = Simulation_SNR_dB_FBMC(i_SNR); 

        Pn                              = SamplingRate/(SubcarrierSpacing*L)*10^(-SNR_dB_FBMC/10)*(K+1)/K; 
        noise_Antenna1                  = sqrt(Pn/2)*( randn(size(s_OFDM_Alamouti_Antenna1)) + 1j*randn(size(s_OFDM_Alamouti_Antenna1)) );
        noise_Antenna2                  = sqrt(Pn/2)*( randn(size(s_OFDM_Alamouti_Antenna1)) + 1j*randn(size(s_OFDM_Alamouti_Antenna1)) );

        %% Received Signal
        r_OFDM_Alamouti_Antenna1        = r_OFDM_Alamouti_noNoise(:,1) + noise_Antenna1; % only Antenna 1 for Alamouti
        r_FBMC_Alamouti_Antenna1        = r_FBMC_Alamouti_noNoise(:,1) + noise_Antenna1; % only Antenna 1 for Alamouti

        r_OFDM_SM_Antenna1              = r_OFDM_SM_noNoise(:,1) + noise_Antenna1;
        r_OFDM_SM_Antenna2              = r_OFDM_SM_noNoise(:,2) + noise_Antenna2;

        r_FBMC_SM_Antenna1              = r_FBMC_SM_noNoise(:,1) + noise_Antenna1;
        r_FBMC_SM_Antenna2              = r_FBMC_SM_noNoise(:,2) + noise_Antenna2;


        %% Demodulation
        y_OFDM_Alamouti_3Blocks         = OFDM.Demodulation(r_OFDM_Alamouti_Antenna1);
        y_FBMC_Alamouti_3Blocks         = FBMC.Demodulation(r_FBMC_Alamouti_Antenna1);
        y_OFDM_Alamouti                 = y_OFDM_Alamouti_3Blocks(:,(1:K/2)+K/2);
        y_FBMC_Alamouti                 = reshape(CBlock'*reshape(y_FBMC_Alamouti_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);

        y_OFDM_SM_Antenna1_3Blocks      = OFDM.Demodulation(r_OFDM_SM_Antenna1);
        y_OFDM_SM_Antenna1              = y_OFDM_SM_Antenna1_3Blocks(:,(1:K/2)+K/2);
        y_OFDM_SM_Antenna2_3Blocks      = OFDM.Demodulation(r_OFDM_SM_Antenna2);
        y_OFDM_SM_Antenna2              = y_OFDM_SM_Antenna2_3Blocks(:,(1:K/2)+K/2);

        y_FBMC_SM_Antenna1_3Blocks      = FBMC.Demodulation(r_FBMC_SM_Antenna1);
        y_FBMC_SM_Antenna1              = reshape(CBlock'*reshape(y_FBMC_SM_Antenna1_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);
        y_FBMC_SM_Antenna2_3Blocks      = FBMC.Demodulation(r_FBMC_SM_Antenna2);
        y_FBMC_SM_Antenna2              = reshape(CBlock'*reshape(y_FBMC_SM_Antenna2_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);


        %% Noise and Power Estimation 
        % Required for real-world measurements but not in simulations
        noise_OFDM_Antenna1_3Blocks     = OFDM.Demodulation( noise_Antenna1 );
        noise_FBMC_Antenna1_3Blocks     = FBMC.Demodulation( noise_Antenna1 );
        noise_OFDM_Antenna1             = noise_OFDM_Antenna1_3Blocks(:, (1:K/2)+K/2 );
        noise_FBMC_Antenna1             = reshape(CBlock' * reshape( noise_FBMC_Antenna1_3Blocks( : , (1:(K+1))+K+1 ) , L*(K+1) , 1 ) , L , K/2 );
        noise_OFDM_Antenna2_3Blocks     = OFDM.Demodulation( noise_Antenna2 );
        noise_FBMC_Antenna2_3Blocks     = FBMC.Demodulation( noise_Antenna2 );
        noise_OFDM_Antenna2             = noise_OFDM_Antenna2_3Blocks(:, (1:K/2 )+ K/2 );
        noise_FBMC_Antenna2             = reshape(CBlock'*reshape(noise_FBMC_Antenna2_3Blocks(:,(1:(K+1))+K+1),L*(K+1),1),L,K/2);

        Pn_OFDM_Antenna1(i_SNR,i_rep)   = mean(abs(noise_OFDM_Antenna1(:)).^2);
        Pn_FBMC_Antenna1(i_SNR,i_rep)   = mean(abs(noise_FBMC_Antenna1(:)).^2);
        Pn_OFDM_Antenna2(i_SNR,i_rep)   = mean(abs(noise_OFDM_Antenna2(:)).^2);
        Pn_FBMC_Antenna2(i_SNR,i_rep)   = mean(abs(noise_FBMC_Antenna2(:)).^2);

        PSignalPlusNoise_OFDM_Alamouti(i_SNR,i_rep)     = mean(abs(y_OFDM_Alamouti(:)).^2);
        PSignalPlusNoise_FBMC_Alamouti(i_SNR,i_rep)     = mean(abs(y_FBMC_Alamouti(:)).^2);

        PSignalPlusNoise_OFDM_SM_Antenna1(i_SNR,i_rep)  = mean(abs(y_OFDM_SM_Antenna1(:)).^2);
        PSignalPlusNoise_OFDM_SM_Antenna2(i_SNR,i_rep)  = mean(abs(y_OFDM_SM_Antenna2(:)).^2);
        PSignalPlusNoise_FBMC_SM_Antenna1(i_SNR,i_rep)  = mean(abs(y_FBMC_SM_Antenna1(:)).^2);
        PSignalPlusNoise_FBMC_SM_Antenna2(i_SNR,i_rep)  = mean(abs(y_FBMC_SM_Antenna2(:)).^2);


        %% Channel Estimation 
        % We average over all pilots to increase the channel estimation accuracy.
        % This only works for a sufficiently low delay spread and a sufficiently 
        % low Doppler spread. In case this is not true, we have to rewrite the
        % channel estimation method. However this averaging is what we used in our
        % 2.5 GHz measurements! 
        EstimatedChannel_OFDM_Alamouti(:,:,1)   = mean(y_OFDM_Alamouti(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1)*ones(L,K/2);
        EstimatedChannel_OFDM_Alamouti(:,:,2)   = mean(y_OFDM_Alamouti(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2)*ones(L,K/2);

        EstimatedChannel_FBMC_Alamouti(:,:,1)   = mean(y_FBMC_Alamouti(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1)*ones(L,K/2);
        EstimatedChannel_FBMC_Alamouti(:,:,2)   = mean(y_FBMC_Alamouti(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2)*ones(L,K/2);

        H11_Est_OFDM    = mean(y_OFDM_SM_Antenna1(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
        H21_Est_OFDM    = mean(y_OFDM_SM_Antenna2(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
        H12_Est_OFDM    = mean(y_OFDM_SM_Antenna1(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);
        H22_Est_OFDM    = mean(y_OFDM_SM_Antenna2(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);

        H11_Est_FBMC    = mean(y_FBMC_SM_Antenna1(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
        H21_Est_FBMC    = mean(y_FBMC_SM_Antenna2(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
        H12_Est_FBMC    = mean(y_FBMC_SM_Antenna1(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);
        H22_Est_FBMC    = mean(y_FBMC_SM_Antenna2(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);

        H_Est_OFDM      = [H11_Est_OFDM , H12_Est_OFDM ; H21_Est_OFDM , H22_Est_OFDM];
        H_Est_FBMC      = [H11_Est_FBMC , H12_Est_FBMC ; H21_Est_FBMC , H22_Est_FBMC];


        %% Data Detection
        % Alamouti
        x_est_OFDM_Alamouti     = Alamouti.Decoder( y_OFDM_Alamouti , EstimatedChannel_OFDM_Alamouti*sqrt(2) );
        x_est_FBMC_Alamouti     = Alamouti.Decoder( y_FBMC_Alamouti , EstimatedChannel_FBMC_Alamouti*sqrt(2) );

        % ML Detection
        y_OFDM_SM_Temp          = repmat(reshape([ y_OFDM_SM_Antenna1(:).' ; y_OFDM_SM_Antenna2(:).' ],2,1,[]) , 1 , QAM.ModulationOrder^2 , 1 );
        [~,indexMin]            = min(sum(abs( y_OFDM_SM_Temp - repmat(H_Est_OFDM*ML_Mapping,1,1,L*K/2) ).^2 , 1 ),[],2);
        x_est_OFDM_SM_ML        = reshape( ML_Mapping(:,indexMin(:)).' , L , K/2 , 2 );

        y_FBMC_SM_Temp          = repmat(reshape([ y_FBMC_SM_Antenna1(:).' ; y_FBMC_SM_Antenna2(:).'],2,1,[]) , 1 , QAM.ModulationOrder^2 , 1);
        [~,indexMin]            = min(sum(abs( y_FBMC_SM_Temp - repmat(H_Est_FBMC*ML_Mapping,1,1,L*K/2) ).^2 , 1 ),[],2);
        x_est_FBMC_SM_ML        = reshape( ML_Mapping(:,indexMin(:)).' , L , K/2 , 2 );


        % Symbols To Bit
        DetectedBitStream_OFDM_Alamouti     = QAM.Symbol2Bit(x_est_OFDM_Alamouti(PilotMatrixBlockAntenna1==0));
        DetectedBitStream_FBMC_Alamouti     = QAM.Symbol2Bit(x_est_FBMC_Alamouti(PilotMatrixBlockAntenna1==0));

        DetectedBitStream_OFDM_SM_ML        = QAM.Symbol2Bit(x_est_OFDM_SM_ML(PilotMatrixBlockAntenna==0));
        DetectedBitStream_FBMC_SM_ML        = QAM.Symbol2Bit(x_est_FBMC_SM_ML(PilotMatrixBlockAntenna==0));

        % Bit Error Ratio
        BER_OFDM_Alamouti(i_SNR,i_rep)      = mean(BinaryDataStream_Alamouti~=DetectedBitStream_OFDM_Alamouti);
        BER_FBMC_Alamouti(i_SNR,i_rep)      = mean(BinaryDataStream_Alamouti~=DetectedBitStream_FBMC_Alamouti);

        BER_OFDM_SM_ML(i_SNR,i_rep)         = mean([BinaryDataStream_SMAntenna1;BinaryDataStream_SMAntenna2]~=DetectedBitStream_OFDM_SM_ML);
        BER_FBMC_SM_ML(i_SNR,i_rep)         = mean([BinaryDataStream_SMAntenna1;BinaryDataStream_SMAntenna2]~=DetectedBitStream_FBMC_SM_ML);


        %% Data Detection: Perfect CSI
        if Simulation_IncludePerfectCSI

            x_est_OFDM_Alamouti_PerfectCSI  = Alamouti.Decoder(y_OFDM_Alamouti,PerfectChannel_OFDM_Alamouti*sqrt(2));
            x_est_FBMC_Alamouti_PerfectCSI  = Alamouti.Decoder(y_FBMC_Alamouti,PerfectChannel_FBMC_Alamouti*sqrt(2));

            DetectedBitStream_OFDM_Alamouti_PerfectCSI  = QAM.Symbol2Bit(x_est_OFDM_Alamouti_PerfectCSI(PilotMatrixBlockAntenna1==0));
            DetectedBitStream_FBMC_Alamouti_PerfectCSI  = QAM.Symbol2Bit(x_est_FBMC_Alamouti_PerfectCSI(PilotMatrixBlockAntenna1==0));

            BER_OFDM_Alamouti_PerfectCSI(i_SNR,i_rep)   = mean(BinaryDataStream_Alamouti~=DetectedBitStream_OFDM_Alamouti_PerfectCSI);
            BER_FBMC_Alamouti_PerfectCSI(i_SNR,i_rep)   = mean(BinaryDataStream_Alamouti~=DetectedBitStream_FBMC_Alamouti_PerfectCSI);


            x_est_OFDM_SM_ML_PerfectCSI = nan(L*K/2,2);
            x_est_FBMC_SM_ML_PerfectCSI = nan(L*K/2,2); 
            for i_lk = 1:NrTransmittedSymbols
                 % OFDM 
                 y_OFDM_rep  = repmat([ y_OFDM_SM_Antenna1(i_lk) ; y_OFDM_SM_Antenna2(i_lk )] , 1 , size(ML_Mapping,2) );
                 H_OFDM_temp = [ H11_OFDM(i_lk) , H12_OFDM(i_lk) ; H21_OFDM(i_lk) , H22_OFDM(i_lk) ];         
                 [~,indexMin] = min(sum( abs(y_OFDM_rep-H_OFDM_temp * ML_Mapping).^2 ,1),[],2);         
                 x_est_OFDM_SM_ML_PerfectCSI(i_lk,:) = ML_Mapping(:,indexMin(:));

                 % FBMC
                 y_FBMC_rep  = repmat([ y_FBMC_SM_Antenna1(i_lk) ; y_FBMC_SM_Antenna2(i_lk )] , 1 , size(ML_Mapping,2) );
                 H_FBMC_temp = [ H11_FBMC(i_lk) , H12_FBMC(i_lk) ; H21_FBMC(i_lk) , H22_FBMC(i_lk) ];         
                 [~,indexMin] = min(sum( abs(y_FBMC_rep-H_FBMC_temp * ML_Mapping).^2 ,1),[],2);         
                 x_est_FBMC_SM_ML_PerfectCSI(i_lk,:) = ML_Mapping(:,indexMin(:));         

            end
            x_est_OFDM_SM_ML_PerfectCSI = reshape(x_est_OFDM_SM_ML_PerfectCSI,L,K/2,2);
            x_est_FBMC_SM_ML_PerfectCSI = reshape(x_est_FBMC_SM_ML_PerfectCSI,L,K/2,2);

            DetectedBitStream_OFDM_SM_ML_PerfectCSI        = QAM.Symbol2Bit(x_est_OFDM_SM_ML_PerfectCSI(PilotMatrixBlockAntenna==0));
            DetectedBitStream_FBMC_SM_ML_PerfectCSI        = QAM.Symbol2Bit(x_est_FBMC_SM_ML_PerfectCSI(PilotMatrixBlockAntenna==0));

            BER_OFDM_SM_ML_PerfectCSI(i_SNR,i_rep)         = mean([BinaryDataStream_SMAntenna1;BinaryDataStream_SMAntenna2]~=DetectedBitStream_OFDM_SM_ML_PerfectCSI);
            BER_FBMC_SM_ML_PerfectCSI(i_SNR,i_rep)         = mean([BinaryDataStream_SMAntenna1;BinaryDataStream_SMAntenna2]~=DetectedBitStream_FBMC_SM_ML_PerfectCSI);

        end
    end
    TimePassed = toc;
    if mod(i_rep,100)==0
        disp(['Realization ' int2str(i_rep) ' of ' int2str(Simulation_MonteCarloRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(Simulation_MonteCarloRepetitions-i_rep)/60) 'minutes']);
    end
end

Estimated_SNR_dB_OFDM = 10*log10((mean(PSignalPlusNoise_OFDM_Alamouti,2)-mean(Pn_OFDM_Antenna1,2))./mean(Pn_OFDM_Antenna1,2));
Estimated_SNR_dB_FBMC = 10*log10((mean(PSignalPlusNoise_FBMC_Alamouti,2)-mean(Pn_FBMC_Antenna1,2))./mean(Pn_FBMC_Antenna1,2));

% To check if both antennas have the same SNR
Estimated_SNR_dB_OFDM_Antenna1 = 10*log10((mean(PSignalPlusNoise_OFDM_SM_Antenna1,2)-mean(Pn_OFDM_Antenna1,2))./mean(Pn_OFDM_Antenna1,2));
Estimated_SNR_dB_OFDM_Antenna2 = 10*log10((mean(PSignalPlusNoise_OFDM_SM_Antenna2,2)-mean(Pn_OFDM_Antenna2,2))./mean(Pn_OFDM_Antenna2,2));
Estimated_SNR_dB_FBMC_Antenna1 = 10*log10((mean(PSignalPlusNoise_FBMC_SM_Antenna1,2)-mean(Pn_FBMC_Antenna1,2))./mean(Pn_FBMC_Antenna1,2));
Estimated_SNR_dB_FBMC_Antenna2 = 10*log10((mean(PSignalPlusNoise_FBMC_SM_Antenna2,2)-mean(Pn_FBMC_Antenna2,2))./mean(Pn_FBMC_Antenna2,2));

% The SNR in OFDM is slightly smaller because we use a zero guard symbol in
% Coded FBMC-OQAM, which requires no power.
M_SNR_dB_OFDM = Simulation_SNR_dB_FBMC - 10*log10((K+1)/K);


%% Calculate BER (mean and confidence intervals)
MeanBER_OFDM_Alamouti               = mean( BER_OFDM_Alamouti ,2);
ConInterval_BER_OFDM_Alamouti       = bootci( 2000 , @(x)(mean(x)) , BER_OFDM_Alamouti' );
ConInterval_L_BER_OFDM_Alamouti     = MeanBER_OFDM_Alamouti - ConInterval_BER_OFDM_Alamouti(1,:)';
ConInterval_U_BER_OFDM_Alamouti     = ConInterval_BER_OFDM_Alamouti(2,:)' - MeanBER_OFDM_Alamouti;

MeanBER_FBMC_Alamouti               = mean( BER_FBMC_Alamouti ,2);
ConInterval_BER_FBMC_Alamouti       = bootci( 2000 , @(x)(mean(x)) , BER_FBMC_Alamouti' );
ConInterval_L_BER_FBMC_Alamouti     = MeanBER_FBMC_Alamouti - ConInterval_BER_FBMC_Alamouti(1,:)';
ConInterval_U_BER_FBMC_Alamouti     = ConInterval_BER_FBMC_Alamouti(2,:)' - MeanBER_FBMC_Alamouti;

MeanBER_OFDM_SM_ML                  = mean( BER_OFDM_SM_ML ,2);
ConInterval_BER_OFDM_SM_ML          = bootci( 2000 , @(x)(mean(x)) , BER_OFDM_SM_ML' );
ConInterval_L_BER_OFDM_SM_ML        = MeanBER_OFDM_SM_ML - ConInterval_BER_OFDM_SM_ML(1,:)';
ConInterval_U_BER_OFDM_SM_ML        = ConInterval_BER_OFDM_SM_ML(2,:)' - MeanBER_OFDM_SM_ML;

MeanBER_FBMC_SM_ML                  = mean( BER_FBMC_SM_ML ,2);
ConInterval_BER_FBMC_SM_ML          = bootci( 2000 , @(x)(mean(x)) , BER_FBMC_SM_ML' );
ConInterval_L_BER_FBMC_SM_ML        = MeanBER_FBMC_SM_ML - ConInterval_BER_FBMC_SM_ML(1,:)';
ConInterval_U_BER_FBMC_SM_ML        = ConInterval_BER_FBMC_SM_ML(2,:)' - MeanBER_FBMC_SM_ML;

%% Plot Results
figure(1);
errorbar( Simulation_SNR_dB_FBMC , MeanBER_FBMC_Alamouti , ConInterval_L_BER_FBMC_Alamouti , ConInterval_U_BER_FBMC_Alamouti , 'LineStyle','none','Color',[0.65,0.65,0.65]);
hold on;
plot(     Simulation_SNR_dB_FBMC , MeanBER_FBMC_Alamouti ,'blue -x');
plot(     Simulation_SNR_dB_FBMC , MeanBER_OFDM_Alamouti ,'red :o');
errorbar( Simulation_SNR_dB_FBMC , MeanBER_FBMC_SM_ML , ConInterval_L_BER_FBMC_SM_ML , ConInterval_U_BER_FBMC_SM_ML , 'LineStyle','none','Color',[0.65,0.65,0.65]);
p2=plot(  Simulation_SNR_dB_FBMC , MeanBER_FBMC_SM_ML ,'blue -x');
p1=plot(  Simulation_SNR_dB_FBMC , MeanBER_OFDM_SM_ML ,':o red');
set(gca,'yScale','log');
xlabel('Estimated Signal-to-Noise Ratio for FBMC [dB]');
ylabel('Bit Error Ratio');
xlim([min(Simulation_SNR_dB_FBMC) max(Simulation_SNR_dB_FBMC)])
l1 = legend([p1 p2],'CP-OFDM','Coded FBMC-OQAM','Location','SouthWest');


if Simulation_IncludePerfectCSI
    figure(2);
    semilogy( Simulation_SNR_dB_FBMC , mean(BER_FBMC_Alamouti_PerfectCSI,2) ,'black -x'); hold on;
    semilogy( Simulation_SNR_dB_FBMC , mean(BER_OFDM_Alamouti_PerfectCSI,2) ,'black :o');
    semilogy( Simulation_SNR_dB_FBMC , MeanBER_FBMC_Alamouti ,'blue -x');
    semilogy( Simulation_SNR_dB_FBMC , MeanBER_OFDM_Alamouti ,'red :o');
    semilogy(  Simulation_SNR_dB_FBMC , mean(BER_FBMC_SM_ML_PerfectCSI,2) ,'black -x');
    semilogy(  Simulation_SNR_dB_FBMC , mean(BER_OFDM_SM_ML_PerfectCSI,2) ,':o black');
    semilogy(  Simulation_SNR_dB_FBMC , MeanBER_FBMC_SM_ML ,'blue -x');
    semilogy(  Simulation_SNR_dB_FBMC , MeanBER_OFDM_SM_ML ,':o red');
    xlabel('Estimated Signal-to-Noise Ratio for FBMC [dB]');
    ylabel('Bit Error Ratio');
    title('Perfect CSI (Black Curve)');
end







