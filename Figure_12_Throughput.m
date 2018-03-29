% =========================================================================   
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =========================================================================      
% This script emulates our throughput measurements. After Turbo decoding,
% all bits must be correctly detected. If one bit is false, the throughput 
% is zero. We maximize the throughput over 15 CQI values, representing 
% perfect feeback. The channel estimation in FBMC is done either by
%       1) Coding method (same as in our measurements), or
%       2) {1,2,3,4}-Auxiliary symbols 
% Note that the current settings of the channel interpolation achieves only
% a good performance for a low delay spread and a low doppler spread. 
% If one wants to test other channel conditions, the interpolation method 
% must be changed, e.g., set to 'linear', or perfect channel knowledge must
% be considered: "Simulation_IncludePerfectCSI = true;".
% For more information about channel estimation in FBMC, we refer to 
% "On Pilot-Symbol Aided Channel Estimation in FBMC-OQAM", R.Nissel, et.al.
% !!!! We recommend using parfor, that is, commenting out line 198-200 !!!!


clear; close all;


%% Parameters
% Simulation parameters
Simulation_MonteCarloRepetitions            = 10;                       % Number of Monte Carlo repetitions (must be larger than 2 so that bootci works)
Simulation_SNR_FBMC_dB                      = -10:5:30;                 % Signal-to-noise ratio for FBMC in dB
Simulation_IncludePerfectCSI                = false;                    % If set to true, also calculate the throughput for perfect channel knowledge (takes more time)
Simulation_PlotAdditionalInformation        = true;                     % If set to true, plot additional information: expected signal power in time, power spectral density, pilot pattern 

% Channel parameters
Channel_PowerDelayProfile                   = 'PedestrianA';            % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)     
Channel_Velocity_kmh                        = 0;                        % Velocity in km/h 
Channel_CarrierFrequency                    = 2.5e9;                    % Carrier frequency (has no influence if the velocity is 0) 

% Modulation parameters
SubcarrierSpacing                           = 15e3;                     % Subcarrier spacing (15kHz, same as LTE)

FBMC_NumberOfSubcarriers                    = 87;                       % The number of subcarriers in FBMC. Can be higher than in OFDM due to lower out of band emissions. 1.4MHz/15kHz=93.33. We consider 6 guard symbols => L=93-6=87. We could also choose L=90 (however, it should be a multiple of 3 so that the pilot density is the same as in OFDM).
FBMC_NumberOfSymbolsInTime                  = 30;                       % Number of FBMC symbols in time. 30*1/15e3/2=1ms!
FBMC_ImaginaryInterferenceCancelationMethod = 'Coding';                 % 'Coding' or 'Auxiliary': Choose between two different methods to cancel the imaginary interference at pilot positions in FBMC
FBMC_NumberOfAuxilarySymbolsPerPilot        = 1;                        % Only relevant if the cancellation method is 'Auxiliary'. Defines the number of auxiliary symbols per pilot. For example, two auxiliary symbols are better than one auxiliary symbol in a low SNR regime. Additionally, the PAPR is better. See "Experimental Evaluation of FBMC-OQAM Channel Estimation Based on Multiple Auxiliary Symbols" R.Nissel et.al

OFDM_NumberOfSubcarriers                    = 72;                       % 1.4 MHz LTE => 72 subcarrier
OFDM_NumberOfSymbolsInTime                  = 14;                       % Number of OFDM symbols (one subframe). LTE uses 14.
OFDM_CyclicPrefixLength                     = 1/(14*15e3);              % LTE CP length: 4.76µs

SamplingRate = SubcarrierSpacing*12*14*4;                               % Multiple of 14 so that CP fits into sampling rate

%% CQI Table
% The first column represents the modulation order: 4, 16, 64, 256, 1024...
% The second column represents the code rate (must be between zero and one)
% Currently, values are chosen according to the (old) LTE standard:
M_CQI = [4  ,  78/1024;...
         4  , 120/1024;...
         4  , 193/1024;...
         4  , 308/1024;...
         4  , 449/1024;...
         4  , 602/1024;...
         16 , 378/1024;...
         16 , 490/1024;...
         16 , 616/1024;...
         64 , 466/1024;...
         64 , 567/1024;...
         64 , 666/1024;...
         64 , 772/1024;...
         64 , 873/1024;...
         64 , 948/1024]; % page 48 of http://www.etsi.org/deliver/etsi_ts/136200_136299/136213/08.08.00_60/ts_136213v080800p.pdf 
     
if not(strcmp(mexext,'mexw64'))  
    % We use a win64 mexfile for code rates smaller than 1/3 => only works 
    % in 64-bit Windows
    IndexCodeRateSmallerOneThird =  find(M_CQI(:,2)<1/3);
    if  numel(IndexCodeRateSmallerOneThird)>0
        M_CQI(IndexCodeRateSmallerOneThird,:) = [];
        warning('A code rate smaller than 1/3 is only supported for Windows 64-bit => CQI values which contain a code rate smaller than 1/3 are discarded!');
    end    
end

     
%% Objects
% FBMC Object
FBMC = Modulation.FBMC(...
    FBMC_NumberOfSubcarriers,...                                        % Number of subcarriers
    FBMC_NumberOfSymbolsInTime,...                                      % Number FBMC symbols in time
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz).  Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    'Hermite-OQAM',...                                                  % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM. The data rate of QAM is reduced by a factor of two compared to OQAM, but robustness in doubly-selective channels is inceased
    4, ...                                                              % Overlapping factor (also determines oversampling in the frequency domain)                                   
    0, ...                                                              % Initial phase shift
    true ...                                                            % Polyphase implementation
    );
FBMC_BlockOverlapTime = (FBMC.PrototypeFilter.OverlappingFactor-1/2)*FBMC.PHY.TimeSpacing;

% OFDM Object
OFDM = Modulation.OFDM(...
    OFDM_NumberOfSubcarriers,...                                        % Number of subcarriers
    OFDM_NumberOfSymbolsInTime,...                                      % Number OFDM symbols in time                                                 
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz) 
    SamplingRate,...                                                    % Sampling rate (Samples/s)                                       
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    OFDM_CyclicPrefixLength, ...                                        % Length of the cyclic prefix (s)                 
    FBMC_BlockOverlapTime ...                                           % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    );

if FBMC.Nr.SamplesTotal == OFDM.Nr.SamplesTotal
    N = FBMC.Nr.SamplesTotal;
else
    error('Number of samples in OFDM and FBMC have to be the same.');   % To simplify the evaluation
end

% Channel object
ChannelModel = Channel.FastFading(...
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    Channel_PowerDelayProfile,...                                       % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate) 
    N,...                                                               % Number of total samples
    Channel_Velocity_kmh/3.6*Channel_CarrierFrequency/2.998e8,...       % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8  
    'Jakes',...                                                         % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large                                       
    200, ...                                                            % Number of paths for the WSSUS process. Only relevant for a 'Jakes' and 'Uniform' Doppler spectrum                                                 
    1,...                                                               % Number of transmit antennas
    1,...                                                               % Number of receive antennas
    true ...                                                            % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );

% Channel Estimation Objects
ChannelEstimation_OFDM = ChannelEstimation.PilotSymbolAidedChannelEstimation(...
    'Diamond',...                                                       % Pilot pattern, 'Diamond','Rectangular', 'Custom'
    [...                                                                % Matrix that represents the pilot pattern parameters
    OFDM.Nr.Subcarriers,...                                                 % Number of subcarriers
    6; ...                                                                  % Pilot spacing in the frequency domain
    OFDM.Nr.MCSymbols,...                                                   % Number of OFDM Symbols
    3.5 ...                                                                 % Pilot spacing in the time domain
    ],...                                   
    'MovingBlockAverage' ,...                                           % Interpolation method: 'MovingBlockAverage' (takes the average of a few close pilots to estimate the channel at the data position), 'FullAverage','linear','nearest','natural'
    [6 OFDM.Nr.MCSymbols] ...                                           % For 'MovingBlockAverage' defines the average region in frequency and time    
    );
ChannelEstimation_FBMC = ChannelEstimation.PilotSymbolAidedChannelEstimation(...
    'Diamond',...                                                       % Pilot pattern, 'Diamond', 'Rectangular', 'Custom'
    [...                                                                % Matrix that represents the pilot pattern parameters
    FBMC.Nr.Subcarriers,...                                                 % Number of subcarriers
    6; ...                                                                  % Pilot spacing in the frequency domain
    FBMC.Nr.MCSymbols,...                                                   % Number of FBMC Symbols
    8 ...                                                                   % Pilot spacing in the time domain
    ],...                                   
    'MovingBlockAverage', ...                                           % Interpolation method: 'MovingBlockAverage' (takes the average of a few close pilots to estimate the channel at the data position), 'FullAverage','linear','nearest','natural'
    [6 FBMC.Nr.MCSymbols]);                                             % For 'MovingBlockAverage' defines the average region in frequency and time    

% FBMC: imaginary interference cancellation method at pilot position
ImagInterCancel  = ChannelEstimation.ImaginaryInterferenceCancellationAtPilotPosition(...
    FBMC_ImaginaryInterferenceCancelationMethod, ...                                        % Cancellation method, either 'Coding' or 'Auxiliary'
    ChannelEstimation_FBMC.GetAuxiliaryMatrix(FBMC_NumberOfAuxilarySymbolsPerPilot), ...    % PilotMatrix, 0 = Data, 1 = Pilot, -1 = Auxiliary symbol  
    FBMC.GetFBMCMatrix, ...                                                                 % FBMC transmission matrix D, i.e., y = D*x with x transmitted data symbols and y received data symbols (before equalization)
    16, ...                                                                                 % Cancel 16 closest interferers
    2 ...                                                                                   % Pilot to data power offset. 2 guarantees that the SNR is the same at pilot position and at data position => fair comparision.
    );


%% Pre-calculate
NrPilotSymbols_OFDM     = ChannelEstimation_OFDM.NrPilotSymbols;
NrDataSymbols_OFDM      = OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols-NrPilotSymbols_OFDM;
NrPilotSymbols_FBMC     = ImagInterCancel.NrPilotSymbols;
NrDataSymbols_FBMC      = ImagInterCancel.NrDataSymbols;

% Pre-initialize CQI: Turbo Coder and QAM
for i_cqi = 1:size(M_CQI,1)
    QAMModulationOrder  = M_CQI(i_cqi,1);
    PAMModulationOrder  = sqrt(QAMModulationOrder);
    CodeRate            = M_CQI(i_cqi,2);
    
    QAM{i_cqi}          = Modulation.SignalConstellation(QAMModulationOrder,'QAM');
    PAM{i_cqi}          = Modulation.SignalConstellation(PAMModulationOrder,'PAM');

    OFDM_TurboCoding{i_cqi} = Coding.TurboCoding(log2(QAMModulationOrder)*NrDataSymbols_OFDM,round(CodeRate*log2(QAMModulationOrder)*NrDataSymbols_OFDM));
    FBMC_TurboCoding{i_cqi} = Coding.TurboCoding(log2(PAMModulationOrder)*NrDataSymbols_FBMC,round(CodeRate*log2(PAMModulationOrder)*NrDataSymbols_FBMC));    
end

% Calculate the SNR. The per-symbol SNR in FBMC and OFDM is different due
% to a different bandwidth and because energy might be wasted (auxiliary
% symbols)
M_SNR_FBMC_dB = Simulation_SNR_FBMC_dB;
for i_SNR = 1:length(Simulation_SNR_FBMC_dB)
    Pn_time                 = 1/FBMC.GetSymbolNoisePower(1) * 10^(-Simulation_SNR_FBMC_dB(i_SNR)/10) * 2 * ImagInterCancel.DataPowerReduction;
    M_SNR_OFDM_dB(i_SNR)    = -10*log10(OFDM.GetSymbolNoisePower(Pn_time));
end

% Preallocate simulation results (needed for parfor)
M_Througput_OFDM            = nan( length(M_SNR_OFDM_dB) , Simulation_MonteCarloRepetitions , size(M_CQI,1) );
M_Througput_OFDM_PerfectCSI = nan( length(M_SNR_OFDM_dB) , Simulation_MonteCarloRepetitions , size(M_CQI,1) );
M_Througput_FBMC            = nan( length(M_SNR_OFDM_dB) , Simulation_MonteCarloRepetitions , size(M_CQI,1) );
M_Througput_FBMC_PerfectCSI = nan( length(M_SNR_OFDM_dB) , Simulation_MonteCarloRepetitions , size(M_CQI,1) );

NrWorkers = 1;                                                          % Conventional FOR loop                   
for i_Rep = 1:Simulation_MonteCarloRepetitions                          % Conventional FOR loop  
% cluster     = parcluster('local');                                    % PARFOR
% NrWorkers   = cluster.NumWorkers;                                     % PARFOR 
% parfor i_Rep = 1:Simulation_MonteCarloRepetitions                     % PARFOR 
    tic
    % Set random channel and noise
    ChannelModel.NewRealization;
    noise_Unscaled = sqrt(1/2)*(randn(N,1)+1j*randn(N,1));

    % Perfect CSI (only an approximation because pulse form is not taken 
    % into account. However, the error is very small)
    h_perfect_OFDM  = ChannelModel.GetTransferFunction( OFDM.GetTimeIndexMidPos, OFDM.Implementation.FFTSize , (1:OFDM.Nr.Subcarriers)+OFDM.Implementation.IntermediateFrequency );
    h_perfect_FBMC  = ChannelModel.GetTransferFunction( FBMC.GetTimeIndexMidPos, FBMC.Implementation.FFTSize , (1:FBMC.Nr.Subcarriers)+FBMC.Implementation.IntermediateFrequency );

    % Preallocate simulation result for one realization (needed for parfor)
    M_Througput_OFDM_OneRealization             = nan( length(M_SNR_OFDM_dB) , size(M_CQI,1) );
    M_Througput_OFDM_PerfectCSI_OneRealization  = nan( length(M_SNR_OFDM_dB) , size(M_CQI,1) );
    M_Througput_FBMC_OneRealization             = nan( length(M_SNR_OFDM_dB) , size(M_CQI,1) );
    M_Througput_FBMC_PerfectCSI_OneRealization  = nan( length(M_SNR_OFDM_dB) , size(M_CQI,1) );

    for i_cqi = 1:size(M_CQI,1)
        % Generate Data stream
        BinaryDataStream_OFDM = randi( [0 1] , OFDM_TurboCoding{i_cqi}.NrDataBits , 1 );
        BinaryDataStream_FBMC = randi( [0 1] , FBMC_TurboCoding{i_cqi}.NrDataBits , 1 );

        % Encode Data stream
        OFDM_TurboCoding{i_cqi}.UpdateInterleaving;
        FBMC_TurboCoding{i_cqi}.UpdateInterleaving;

        CodedBits_OFDM  = OFDM_TurboCoding{i_cqi}.TurboEncoder( BinaryDataStream_OFDM );
        CodedBits_FBMC  = FBMC_TurboCoding{i_cqi}.TurboEncoder( BinaryDataStream_FBMC );
        
        % Bit Interleaving
        BitInterleaving_OFDM = randperm( OFDM_TurboCoding{i_cqi}.NrCodedBits );
        BitInterleaving_FBMC = randperm( FBMC_TurboCoding{i_cqi}.NrCodedBits );

        CodedBits_OFDM = CodedBits_OFDM( BitInterleaving_OFDM );
        CodedBits_FBMC = CodedBits_FBMC( BitInterleaving_FBMC );

        % Pilot Symbols
        xP_OFDM = QAM{i_cqi}.SymbolMapping( randi(QAM{i_cqi}.ModulationOrder,[NrPilotSymbols_OFDM 1]) );
        xP_OFDM = xP_OFDM./abs(xP_OFDM);
        xP_FBMC = PAM{i_cqi}.SymbolMapping( randi(PAM{i_cqi}.ModulationOrder,[NrPilotSymbols_FBMC 1]) );
        xP_FBMC = xP_FBMC./abs(xP_FBMC);

        % Map coded bits to symbols
        xD_OFDM = QAM{i_cqi}.Bit2Symbol( CodedBits_OFDM );
        xD_FBMC = PAM{i_cqi}.Bit2Symbol( CodedBits_FBMC );

        % transmitted symbols
        x_OFDM                                          = nan(OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols);
        x_OFDM(ChannelEstimation_OFDM.PilotMatrix==1)   = xP_OFDM;
        x_OFDM(ChannelEstimation_OFDM.PilotMatrix==0)   = xD_OFDM;

        x_FBMC = reshape(ImagInterCancel.PrecodingMatrix*[xP_FBMC;xD_FBMC],[FBMC.Nr.Subcarriers FBMC.Nr.MCSymbols]);

        % Generate the transmitted OFDM and FBMC signal in time
        s_OFDM  = OFDM.Modulation(x_OFDM);
        s_FBMC  = FBMC.Modulation(x_FBMC);

        % Received signal without noise
        r_OFDM_noNoise  = ChannelModel.Convolution(s_OFDM);
        r_FBMC_noNoise  = ChannelModel.Convolution(s_FBMC);

        for i_SNR = 1:length(M_SNR_OFDM_dB)
            SNR_OFDM_dB = M_SNR_OFDM_dB(i_SNR);
            Pn_time     = 1/OFDM.GetSymbolNoisePower(1)*10^(-SNR_OFDM_dB/10);

            % For unit data symbol power and constant transmit power, the noise power scales accordingly
            Pn_OFDM = OFDM.GetSymbolNoisePower( Pn_time );
            Pn_FBMC = FBMC.GetSymbolNoisePower( Pn_time );

            % Add white Gaussian noise
            noise   = sqrt(Pn_time)*noise_Unscaled;
            r_OFDM  = r_OFDM_noNoise + noise;
            r_FBMC  = r_FBMC_noNoise + noise;

            % Demodulate OFDM and FBMC signal
            y_OFDM  = OFDM.Demodulation(r_OFDM);
            y_FBMC  = FBMC.Demodulation(r_FBMC);

            % LS channel estimates at the pilot positions
            hP_est_OFDM = y_OFDM(ChannelEstimation_OFDM.PilotMatrix==1)./xP_OFDM;
            hP_est_FBMC = y_FBMC(ChannelEstimation_FBMC.PilotMatrix==1)./xP_FBMC/sqrt(ImagInterCancel.PilotToDataPowerOffset*ImagInterCancel.DataPowerReduction);

            % Channel estimation (interpolation)
            h_est_OFDM  = reshape( ChannelEstimation_OFDM.ChannelInterpolation(hP_est_OFDM) , [OFDM.Nr.Subcarriers OFDM.Nr.MCSymbols] );
            h_est_FBMC  = reshape( ChannelEstimation_FBMC.ChannelInterpolation(hP_est_FBMC) , [FBMC.Nr.Subcarriers FBMC.Nr.MCSymbols] );

            % Equalize received symbols at data position
            y_EQ_OFDM                           = y_OFDM(ChannelEstimation_OFDM.PilotMatrix==0) ./ h_est_OFDM( ChannelEstimation_OFDM.PilotMatrix==0  );
            NoiseScaling_OFDM                   = 1./abs(h_est_OFDM(ChannelEstimation_OFDM.PilotMatrix==0 )).^2;
                       
            y_EQ_FBMC_all                       = y_FBMC ./ h_est_FBMC;
            if strcmp(FBMC_ImaginaryInterferenceCancelationMethod,'Coding')
                % Coding method: we need despreading at the receiver
                y_EQ_FBMC                       = ImagInterCancel.PrecodingMatrix( : , ImagInterCancel.NrPilotSymbols+1:end )' * y_EQ_FBMC_all(:);
                h_temp                          = ImagInterCancel.PostCodingChannelMatrix(ImagInterCancel.NrPilotSymbols+1:end,:)*h_est_FBMC(:);
                NoiseScaling_FBMC               = 1./abs(h_temp).^2;
            else
                % Auxiliary symbols method
                y_EQ_FBMC                       = y_EQ_FBMC_all( ImagInterCancel.PilotMatrix==0 ) / sqrt( ImagInterCancel.DataPowerReduction );    
                NoiseScaling_FBMC               = 1./abs(h_est_FBMC(  ImagInterCancel.PilotMatrix==0 )).^2;    
            end            

            % Calculate LLR assuming perfect channel knowledge and ignoring interference
            LLR_OFDM    = QAM{i_cqi}.LLR_AWGN( y_EQ_OFDM       , Pn_OFDM .* NoiseScaling_OFDM);
            LLR_FBMC    = PAM{i_cqi}.LLR_AWGN( real(y_EQ_FBMC) , Pn_FBMC .* NoiseScaling_FBMC);

            % Bitdeinterleaving
            LLR_OFDM(BitInterleaving_OFDM)  = LLR_OFDM;
            LLR_FBMC(BitInterleaving_FBMC)  = LLR_FBMC;

            % Decode Bits
            DecodedBits_OFDM    = OFDM_TurboCoding{i_cqi}.TurboDecoder( LLR_OFDM );
            DecodedBits_FBMC    = FBMC_TurboCoding{i_cqi}.TurboDecoder( LLR_FBMC );
                          
            % Simulated throughput after decoding (all bits must be correctly detected. If one bit is wrong, the throughput is zero)
            M_Througput_OFDM_OneRealization(i_SNR,i_cqi)  = all( DecodedBits_OFDM  == BinaryDataStream_OFDM ) * length(BinaryDataStream_OFDM)/(OFDM.PHY.TimeSpacing*(OFDM.Nr.MCSymbols));
            M_Througput_FBMC_OneRealization(i_SNR,i_cqi)  = all( DecodedBits_FBMC  == BinaryDataStream_FBMC ) * length(BinaryDataStream_FBMC)/(FBMC.PHY.TimeSpacing*(FBMC.Nr.MCSymbols));
            
            % Calculate the throughput for perfect channel state information
            if Simulation_IncludePerfectCSI             
                y_EQ_OFDM_PerfectCSI            = y_OFDM(ChannelEstimation_OFDM.PilotMatrix==0) ./ h_perfect_OFDM( ChannelEstimation_OFDM.PilotMatrix==0 );
                NoiseScaling_OFDM_PerfectCSI    = 1./abs(h_perfect_OFDM(ChannelEstimation_OFDM.PilotMatrix==0 )).^2;

                y_EQ_FBMC_all_PerfectCSI    = y_FBMC ./ h_perfect_FBMC;
                if strcmp(FBMC_ImaginaryInterferenceCancelationMethod, 'Coding' )
                    % Coding method: we need despreading at the receiver
                    y_EQ_FBMC_PerfectCSI            = ImagInterCancel.PrecodingMatrix( : , ImagInterCancel.NrPilotSymbols+1:end )' * y_EQ_FBMC_all_PerfectCSI(:);
                    h_temp_PerfectCSI               = ImagInterCancel.PostCodingChannelMatrix(ImagInterCancel.NrPilotSymbols+1:end,:)*h_perfect_FBMC(:);
                    NoiseScaling_FBMC_PerfectCSI    = 1./abs(h_temp_PerfectCSI).^2;
                else
                    % Auxiliary symbols
                    y_EQ_FBMC_PerfectCSI            = y_EQ_FBMC_all_PerfectCSI( ImagInterCancel.PilotMatrix==0 ) / sqrt( ImagInterCancel.DataPowerReduction );    
                    NoiseScaling_FBMC_PerfectCSI    = 1./abs(h_perfect_FBMC(  ImagInterCancel.PilotMatrix==0 )).^2; 
                end

                LLR_OFDM_PerfectCSI     = QAM{i_cqi}.LLR_AWGN( y_EQ_OFDM_PerfectCSI        , Pn_OFDM .* NoiseScaling_OFDM_PerfectCSI);
                LLR_FBMC_PerfectCSI     = PAM{i_cqi}.LLR_AWGN( real(y_EQ_FBMC_PerfectCSI)  , Pn_FBMC .* NoiseScaling_FBMC_PerfectCSI);

                LLR_OFDM_PerfectCSI(BitInterleaving_OFDM)   = LLR_OFDM_PerfectCSI;
                LLR_FBMC_PerfectCSI(BitInterleaving_FBMC)   = LLR_FBMC_PerfectCSI;

                DecodedBits_OFDM_PerfectCSI     = OFDM_TurboCoding{i_cqi}.TurboDecoder(LLR_OFDM_PerfectCSI);
                DecodedBits_FBMC_PerfectCSI     = FBMC_TurboCoding{i_cqi}.TurboDecoder(LLR_FBMC_PerfectCSI);

                M_Througput_OFDM_PerfectCSI_OneRealization(i_SNR,i_cqi) = all( DecodedBits_OFDM_PerfectCSI == BinaryDataStream_OFDM ) * length(BinaryDataStream_OFDM)/(OFDM.PHY.TimeSpacing*(OFDM.Nr.MCSymbols));
                M_Througput_FBMC_PerfectCSI_OneRealization(i_SNR,i_cqi) = all( DecodedBits_FBMC_PerfectCSI == BinaryDataStream_FBMC ) * length(BinaryDataStream_FBMC)/(FBMC.PHY.TimeSpacing*(FBMC.Nr.MCSymbols));
            end
        end
    end

    M_Througput_OFDM(:,i_Rep,:)             = M_Througput_OFDM_OneRealization;
    M_Througput_FBMC(:,i_Rep,:)             = M_Througput_FBMC_OneRealization;
    
    if Simulation_IncludePerfectCSI          
        M_Througput_OFDM_PerfectCSI(:,i_Rep,:)  = M_Througput_OFDM_PerfectCSI_OneRealization;
        M_Througput_FBMC_PerfectCSI(:,i_Rep,:)  = M_Througput_FBMC_PerfectCSI_OneRealization;
    end

    TimePassed = toc;
    disp(['Realization ' int2str(i_Rep) ' of ' int2str(Simulation_MonteCarloRepetitions) 'needed ' int2str(TimePassed) 's. Total simulation time:' int2str(TimePassed*Simulation_MonteCarloRepetitions/NrWorkers/60) 'minutes']);
end

%% Calculate Maximum Throughput 
% Here, we maximize the throughput over all CQI values, that is, we assume 
% perfect feedback. 
Througput_OFDM                              = max(M_Througput_OFDM,[],3);
Mean_Througput_OFDM                         = mean(Througput_OFDM,2);
ConInterval_Througput_OFDM                  = bootci(2000,@(x)(mean(x)),Througput_OFDM');
ConInterval_L_Througput_OFDM                = Mean_Througput_OFDM - ConInterval_Througput_OFDM(1,:)';
ConInterval_U_Througput_OFDM                = ConInterval_Througput_OFDM(2,:)'-Mean_Througput_OFDM;

Througput_FBMC                              = max(M_Througput_FBMC,[],3);
Mean_Througput_FBMC                         = mean(Througput_FBMC,2);
ConInterval_Througput_FBMC                  = bootci(2000,@(x)(mean(x)),Througput_FBMC');
ConInterval_L_Througput_FBMC                = Mean_Througput_FBMC - ConInterval_Througput_FBMC(1,:)';
ConInterval_U_Througput_FBMC                = ConInterval_Througput_FBMC(2,:)'-Mean_Througput_FBMC;

%% Information Theory
% Calculate the achievable rate (Gaussian inputs)
M_SNR_FBMC_dB_MorePoints = linspace(min(M_SNR_FBMC_dB),max(M_SNR_FBMC_dB),100);
for i_SNR = 1:length(M_SNR_FBMC_dB_MorePoints)
    SNR             = 10.^(M_SNR_FBMC_dB_MorePoints(i_SNR)/10);
    fun             = @(h) log2(1+SNR*abs(h).^2).*2.*h.*exp(-(h.^2));
    C_OneSymbol     = integral(fun,0,inf);
    AchievableRate_FBMC(i_SNR)   = NrDataSymbols_FBMC*1/2*C_OneSymbol/(FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols);
end

% Calculate the achievable rate (BICM: 4-OQAM, 16-OQAM,...)
% The BICM capacity per symbol is precalulated. Use the script "./Theory/BICM_Capacity_Rayleigh" to calculate those values
BICM_Capacity_Rayleigh_2_4_8_PAM        = load('./Theory/BICM_Capacity_Rayleigh_2_4_8_PAM.mat');
AchievableRate_FBMC_BICM_2_4_8_PAM      = NrDataSymbols_FBMC * interp1( BICM_Capacity_Rayleigh_2_4_8_PAM.SNR_dB    , BICM_Capacity_Rayleigh_2_4_8_PAM.C_max,M_SNR_FBMC_dB_MorePoints    , 'spline' )/( FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols );

BICM_Capacity_Rayleigh_2_4_8_16_PAM     = load('./Theory/BICM_Capacity_Rayleigh_2_4_8_16_PAM.mat');
AchievableRate_FBMC_BICM_2_4_8_16_PAM   = NrDataSymbols_FBMC * interp1( BICM_Capacity_Rayleigh_2_4_8_16_PAM.SNR_dB , BICM_Capacity_Rayleigh_2_4_8_16_PAM.C_max,M_SNR_FBMC_dB_MorePoints , 'spline' )/( FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols );

%% Plot Results
fprintf('====================================================================================================\n');
fprintf('================================ Basic Settings (without guard)  ===================================\n');
fprintf('               |(complex)TF-Spacing| Bandwidth(FL)|    Time(KT) | SNR rel. to OFDM | Pilot density | \n');
fprintf('OFDM (with CP) |%17.2f  |%8.2f MHz  |%8.2f ms  |%12.2f dB   |    %8.3f   | \n', OFDM.PHY.TimeSpacing*OFDM.PHY.SubcarrierSpacing   , OFDM.PHY.SubcarrierSpacing*OFDM.Nr.Subcarriers/1e6   , OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols/1e-3   ,                  0                , NrPilotSymbols_OFDM/(OFDM.Nr.MCSymbols*OFDM.PHY.TimeSpacing*OFDM.Nr.Subcarriers*OFDM.PHY.SubcarrierSpacing) );
fprintf('FBMC-OQAM      |%17.2f  |%8.2f MHz  |%8.2f ms  |%12.2f dB   |    %8.3f   | \n', FBMC.PHY.TimeSpacing*FBMC.PHY.SubcarrierSpacing*2 , FBMC.PHY.SubcarrierSpacing*FBMC.Nr.Subcarriers/1e6   , FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols/1e-3   , M_SNR_FBMC_dB(1)-M_SNR_OFDM_dB(1) , NrPilotSymbols_FBMC/(FBMC.Nr.MCSymbols*FBMC.PHY.TimeSpacing*FBMC.Nr.Subcarriers*FBMC.PHY.SubcarrierSpacing) );
fprintf('====================================================================================================\n');
fprintf('====================================================================================================\n');


figure(1);
errorbar( M_SNR_FBMC_dB            , Mean_Througput_FBMC/1e6                    , ConInterval_L_Througput_FBMC/1e6 , ConInterval_U_Througput_FBMC/1e6 , 'blue');
hold on;
errorbar( M_SNR_FBMC_dB            , Mean_Througput_OFDM/1e6                    , ConInterval_L_Througput_OFDM/1e6 , ConInterval_U_Througput_OFDM/1e6 , 'red');
semilogy( M_SNR_FBMC_dB_MorePoints , AchievableRate_FBMC/1e6                    , 'Color' , [0,0.4,0.6]);
semilogy( M_SNR_FBMC_dB_MorePoints , AchievableRate_FBMC_BICM_2_4_8_PAM/1e6     , 'Color' , [0,0.2,0.8]);
semilogy( M_SNR_FBMC_dB_MorePoints , AchievableRate_FBMC_BICM_2_4_8_16_PAM/1e6  , '--','Color' , [0,0.2,0.8]);
xlim([min(M_SNR_FBMC_dB) max(M_SNR_FBMC_dB)]);
xlabel('Signal-to-Noise Ratio for FBMC [dB]');
ylabel('Throughput, Achievable Rate [Mbit/s]');


if Simulation_IncludePerfectCSI        
    Througput_OFDM_PerfectCSI                   = max(M_Througput_OFDM_PerfectCSI,[],3);
    Mean_Througput_OFDM_PerfectCSI              = mean(Througput_OFDM_PerfectCSI,2);

    Througput_FBMC_PerfectCSI                   = max(M_Througput_FBMC_PerfectCSI,[],3);
    Mean_Througput_FBMC_PerfectCSI              = mean(Througput_FBMC_PerfectCSI,2);

    figure(2);
    errorbar( M_SNR_FBMC_dB            , Mean_Througput_FBMC/1e6                    , ConInterval_L_Througput_FBMC/1e6 , ConInterval_U_Througput_FBMC/1e6 , 'blue');
    hold on;
    errorbar( M_SNR_FBMC_dB            , Mean_Througput_OFDM/1e6                    , ConInterval_L_Througput_OFDM/1e6 , ConInterval_U_Througput_OFDM/1e6 , 'red');
    semilogy( M_SNR_FBMC_dB_MorePoints , AchievableRate_FBMC/1e6                    , 'Color' , [0,0.4,0.6]);
    semilogy( M_SNR_FBMC_dB_MorePoints , AchievableRate_FBMC_BICM_2_4_8_PAM/1e6     , 'Color' , [0,0.2,0.8]);
    semilogy( M_SNR_FBMC_dB_MorePoints , AchievableRate_FBMC_BICM_2_4_8_16_PAM/1e6  , '--','Color' , [0,0.2,0.8]);
    semilogy( M_SNR_FBMC_dB            , Mean_Througput_OFDM_PerfectCSI/1e6         , 'black');
    semilogy( M_SNR_FBMC_dB            , Mean_Througput_FBMC_PerfectCSI/1e6         , 'black');
    xlim([min(M_SNR_FBMC_dB) max(M_SNR_FBMC_dB)]);
    xlabel('Signal-to-Noise Ratio for FBMC [dB]');
    ylabel('Throughput, Achievable Rate [Mbit/s]');
end


if Simulation_PlotAdditionalInformation 
    % Expected signal power in time (per subframe)
    disp('Calculate expected signal power in time ...');
    [PS_OFDM , t]   = OFDM.PlotTransmitPower;
    PS_FBMC         = FBMC.PlotTransmitPower(ImagInterCancel.PrecodingMatrix*ImagInterCancel.PrecodingMatrix');
    
    figure(3);
    disp('Calculate power spectral density ...');
    plot((t-FBMC_BlockOverlapTime)/1e-3 , PS_OFDM , 'red'); hold on;
    plot((t-FBMC_BlockOverlapTime)/1e-3 , PS_FBMC , 'blue'); 
    xlabel('Time [ms]');
    ylabel('Signal Power');
    
    % Power spectral density (normalization so that Parseval's theorem is fullfilled)
    [PSD_OFDM , f]  = OFDM.PlotPowerSpectralDensity;
    PSD_FBMC        = FBMC.PlotPowerSpectralDensity(ImagInterCancel.PrecodingMatrix*ImagInterCancel.PrecodingMatrix');    
    
    figure(4);
    plot(f/1e6 , 10*log10(PSD_OFDM) ,'red'); hold on;
    plot(f/1e6 , 10*log10(PSD_FBMC) ,'blue'); 
    xlabel('Frequency [MHz]');
    ylabel('Power Spectral Density [dB]');  
    ylim( [-100 2] + 10*log10( max([PSD_OFDM;PSD_FBMC]) ) );
        
    % Pilot pattern for OFDM and FBMC
    figure(5);
    ChannelEstimation_OFDM.PlotPilotPattern;
    title('Pilot Pattern: OFDM');
    figure(6);
    if strcmp(FBMC_ImaginaryInterferenceCancelationMethod,'Coding');
        ChannelEstimation_FBMC.PlotPilotPattern(-(ImagInterCancel.ConsideredInterferenceMatrix<0)+(ImagInterCancel.ConsideredInterferenceMatrix>0))
    else
        ChannelEstimation_FBMC.PlotPilotPattern(ImagInterCancel.PilotMatrix);        
    end
    title('Pilot Pattern: FBMC');
    
end