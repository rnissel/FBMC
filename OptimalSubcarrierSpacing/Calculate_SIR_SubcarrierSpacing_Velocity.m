% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script calculates the Signal-to-Interference Ratio (SIR) for
% different subcarrier spacings and velocities. The results are saved in
% "./Results" and can then be plotted with the script 
% "Plot_SIR_OptimalSubcarrierSpacing", which additionally calculates the
% optimum subcarrier spacing (exhaustive search)

clear; close all;

cd('..');

%% Parameters
Scenario  = 1;                      % Integer number which represents different pre-defined scenarios (channel model, modulation scheme,...), see below
SaveStuff = false;                  % If true, saves the results in "./Results/xxxx", which can later be used to plot the optimal subcarrier spacing

M_Velocity_kmh = [0 250 500];       % Velocies for which the SIR is calculated

% Velocities used in the paper:
% M_Velocity_kmh = [0 2:2:8 10:10:500];
% M_Velocity_kmh = [2:2:8 12 14 16];     

%% Start calculation
for i_Velocity = 1:length(M_Velocity_kmh)
Velocity_kmh = M_Velocity_kmh(i_Velocity);

switch Scenario
    case 1 % Plot 1 Vehicular A
         Method                 = 'OFDM';
         SamplingRate           = 15e3*14*14;                 %Determined by channen model (fits approximatly the channel delay taps)
         PowerDelayProfile      = 'VehicularA';
         M_SubcarrierSpacing    = 1e3*[5:1:60];        
         PrototypeFilter        = '';
         CarrierFrequency       = 2.5e9;
    case 2 % Plot 1 Vehicular A
         Method                 = 'FBMC-OQAM';
         SamplingRate           = 15e3*14*14;               
         PowerDelayProfile      = 'VehicularA';
         M_SubcarrierSpacing    = 1e3*[5:1:60]; 
         PrototypeFilter        = 'Hermite';  
         CarrierFrequency       = 2.5e9;
    case 3 % Plot 1 Vehicular A
         Method                 = 'FBMC-OQAM';
         SamplingRate           = 15e3*14*14;                
         PowerDelayProfile      = 'VehicularA';
         M_SubcarrierSpacing    = 1e3*[5:1:60]; 
         PrototypeFilter        = 'PHYDYAS';
         CarrierFrequency       = 2.5e9;
    case 4 % Pedestrian A
         Method = 'OFDM';
         SamplingRate           = 15e3*14*12*4; 
         PowerDelayProfile      = 'PedestrianA';
         M_SubcarrierSpacing    = 1e3*[5:1:100];
         PrototypeFilter        = '';
         CarrierFrequency       = 2.5e9;         
    case 5 % Pedestrian A
         Method                 = 'FBMC-OQAM';
         SamplingRate           = 15e3*14*12*4;
         PowerDelayProfile      = 'PedestrianA';
         M_SubcarrierSpacing    = 1e3*[15:1:170];  
         PrototypeFilter        = 'Hermite';
         CarrierFrequency       = 2.5e9;
    case 6 % Pedestrian A
         Method                 = 'FBMC-OQAM';
         SamplingRate           = 15e3*14*12*4;
         PowerDelayProfile      = 'PedestrianA';
         M_SubcarrierSpacing    = 1e3*[15:1:170];   
         PrototypeFilter        = 'PHYDYAS';  
         CarrierFrequency       = 2.5e9;
    case 7 % Vehicular B
         Method                 = 'CP-OFDM';
         SamplingRate           = 15e3*14*7; 
         PowerDelayProfile      = 'VehicularB';
         M_SubcarrierSpacing    = 1e3*[5:1:70]; 
         PrototypeFilter        = ''; 
         CarrierFrequency       = 2.5e9;
    case 8 % Vehicular B
         Method                 = 'FBMC-QAM';
         SamplingRate           = 15e3*14*7;
         PowerDelayProfile      = 'VehicularB';
         M_SubcarrierSpacing    = 1e3*[5:1:70]; 
         PrototypeFilter        = 'Hermite';
         CarrierFrequency       = 2.5e9;
    case 9 % Vehicular B
         Method = 'FBMC-QAM';
         SamplingRate           = 15e3*14*7; 
         PowerDelayProfile      = 'VehicularB';
         M_SubcarrierSpacing    = 1e3*[5:1:70];
         PrototypeFilter        = 'Hermite'; 
         CarrierFrequency       = 2.5e9;
    case 10 % Vehicular B
         Method                 = 'FBMC-OQAM';
         SamplingRate           = 15e3*14*7;
         PowerDelayProfile      = 'VehicularB';
         M_SubcarrierSpacing    = 1e3*[5:1:70];   
         PrototypeFilter        = 'Hermite'; 
         CarrierFrequency       = 2.5e9;
   case 11 % Vehicular B
         Method                 = 'OFDM';
         SamplingRate           = 15e3*14*7;
         PowerDelayProfile      = 'VehicularB';
         M_SubcarrierSpacing    = 1e3*[5:1:70]; 
         PrototypeFilter        = '';
         CarrierFrequency       = 2.5e9;
   case 12 % Vehicular B
         Method                 = 'CP-OFDM-LTE'; 
         SamplingRate           = 15e3*14*7;
         PowerDelayProfile      = 'VehicularB';
         M_SubcarrierSpacing    = 1e3*[5:1:30]; 
         PrototypeFilter        = '';
         CarrierFrequency       = 2.5e9;
   case 13 
         Method                 = 'CP-OFDM-007LTE';
         SamplingRate           = 15e3*14*12*4*5;
         PowerDelayProfile      = 'PedestrianA';
         M_SubcarrierSpacing    = 1e3*[170:1:180]; 
         PrototypeFilter        = ''; 
         CarrierFrequency       = 2.5e9;
   case 14 
         Method                 = 'CP-OFDM-007LTE'; 
         SamplingRate           = 15e3*14*14*4;
         PowerDelayProfile      = 'VehicularA';
         M_SubcarrierSpacing    = 1e3*[26:1:50];  
         PrototypeFilter        = ''; 
         CarrierFrequency       = 2.5e9;
   case 15 
         Method                 = 'CP-OFDM-010LTE'; % CP OFDM with TF=1.10!
         SamplingRate           = 15e3*14*7*1;
         PowerDelayProfile      = 'VehicularB';
         M_SubcarrierSpacing    = 1e3*[5:1:20]; 
         PrototypeFilter        = ''; 
         CarrierFrequency       = 2.5e9;
    
    %% New 38.900 Channel Model!
    case 16
         Method                 = 'FBMC-OQAM';
         SamplingRate           = 1/25e-9;
         PowerDelayProfile      = 'TDL-A_30ns';
         M_SubcarrierSpacing    = 1e3*[100:5:1000];  
         PrototypeFilter        = 'Hermite'; 
         CarrierFrequency       = 60e9;      
    case 17
         Method                 = 'FBMC-OQAM';
         SamplingRate           = 1/25e-9;
         PowerDelayProfile      = 'TDL-A_30ns';
         M_SubcarrierSpacing    = 1e3*[100:5:1000];  
         PrototypeFilter        = 'PHYDYAS'; 
         CarrierFrequency       = 60e9;    
    case 18
         Method                 = 'OFDM';
         SamplingRate           = 1/25e-9;
         PowerDelayProfile      = 'TDL-A_30ns';
         M_SubcarrierSpacing    = 1e3*[100:5:1000];  
         PrototypeFilter        = ''; 
         CarrierFrequency       = 60e9;     
    case 19
         Method                 = 'CP-OFDM-007LTE';
         SamplingRate           = 1/25e-9;
         PowerDelayProfile      = 'TDL-A_30ns';
         M_SubcarrierSpacing    = 1e3*[100:5:1000];  
         PrototypeFilter        = ''; 
         CarrierFrequency       = 60e9;  
    case 20
         Method                 = 'FBMC-QAM';
         SamplingRate           = 1/100e-9;
         PowerDelayProfile      = 'TDL-B_900ns';
         M_SubcarrierSpacing    = 1e3*[100:2.5:500];  
         PrototypeFilter        = 'Hermite'; 
         CarrierFrequency       = 60e9;  
    case 21
         Method                 = 'FBMC-QAM';
         SamplingRate           = 1/100e-9;
         PowerDelayProfile      = 'TDL-B_900ns';
         M_SubcarrierSpacing    = 1e3*[100:2.5:500];  
         PrototypeFilter        = 'PHYDYAS'; 
         CarrierFrequency       = 60e9;    
    case 22
         Method                 = 'CP-OFDM';
         SamplingRate           = 1/100e-9;
         PowerDelayProfile      = 'TDL-B_900ns';
         M_SubcarrierSpacing    = 1e3*[100:2.5:500];  
         PrototypeFilter        = ''; 
         CarrierFrequency       = 60e9;
    case 23
         Method                 = 'FBMC-OQAM';
         SamplingRate           = 1/100e-9;
         PowerDelayProfile      = 'TDL-B_900ns';
         M_SubcarrierSpacing    = 1e3*[100:2.5:500];  
         PrototypeFilter        = 'Hermite'; 
         CarrierFrequency       = 60e9;   
    % Windowed and Filtered OFDM. Note that in Filtered OFDM, the SIR depends on the subcarrier index (and not just the 3dB at the edges)
    case 24
         Method                 = 'WOLA';
         SamplingRate           = 15e3*14*14; 
         PowerDelayProfile      = 'VehicularA';
         M_SubcarrierSpacing    = 1e3*[5:1:60]; 
         PrototypeFilter        = '';
         CarrierFrequency       = 2.5e9;     
    case 25
         Method                 = 'FOFDM';
         SamplingRate           = 15e3*14*14; 
         PowerDelayProfile      = 'VehicularA';
         M_SubcarrierSpacing    = 1e3*[5:1:60]; 
         PrototypeFilter        = '';
         CarrierFrequency       = 2.5e9;             
end

% M_RealSubcarrierSpacing = nan(length(M_SubcarrierSpacing));
% SIR_Theory  = nan(length(M_SubcarrierSpacing));
% parfor i_subcarrier = 1:length(M_SubcarrierSpacing)
for i_subcarrier = 1:length(M_SubcarrierSpacing)
    SubcarrierSpacing = M_SubcarrierSpacing(i_subcarrier);
    switch Method
        case 'OFDM'
            Multicarrier = Modulation.OFDM(...
                floor(SamplingRate/SubcarrierSpacing),...       % Number Subcarriers, critically sampled => represents infinitly many subcarriers
                3,...                                           % Number OFDM Symbols
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmitreal valued signal
                0, ...                                          % Cyclic prefix length (s)
                0 ...                                           % Zero guard length (s)
                );
         case 'CP-OFDM'
            Multicarrier = Modulation.OFDM(...
                floor(SamplingRate/SubcarrierSpacing),...       % Number Subcarriers, critically sampled => represents infinitly many subcarriers
                3,...                                           % Number OFDM Symbols
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmitreal valued signal
                1/SubcarrierSpacing, ...                        % Cyclic prefix length (s)
                0 ...                                           % Zero guard length (s)
                );
         case 'CP-OFDM-LTE'
            Multicarrier = Modulation.OFDM(...
                floor(SamplingRate/SubcarrierSpacing),...       % Number Subcarriers, critically sampled => represents infinitly many subcarriers
                3,...                                           % Number OFDM Symbols
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmitreal valued signal
                0.25/SubcarrierSpacing, ...                     % Cyclic prefix length (s). 0.25=Extended CP
                0 ...                                           % Zero guard length (s)
                ); 
         case 'CP-OFDM-007LTE'
            Multicarrier = Modulation.OFDM(...
                floor(SamplingRate/SubcarrierSpacing),...       % Number Subcarriers, critically sampled => represents infinitly many subcarriers
                3,...                                           % Number OFDM Symbols
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmitreal valued signal
                1/14/SubcarrierSpacing, ...                     % Cyclic prefix length (s) LTE: 1/15e3/14;
                0 ...                                           % Zero guard length (s)
                );     
          case 'CP-OFDM-010LTE'
            Multicarrier = Modulation.OFDM(...
                floor(SamplingRate/SubcarrierSpacing),...       % Number Subcarriers, critically sampled => represents infinitly many subcarriers
                3,...                                           % Number OFDM Symbols
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmitreal valued signal
                1/10/SubcarrierSpacing, ...                     % Cyclic prefix length (s) LTE: 1/15e3/14;
                0 ...                                           % Zero guard length (s)
                );                 
         case 'FBMC-OQAM'
            Multicarrier = Modulation.FBMC(...
                7,...                                           % Number subcarriers. Localization => 7 is enough
                7,...                                           % Number FBMC symbols. Localization => 7 is enough
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmit real valued signal
                [PrototypeFilter '-OQAM'],...                   % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
                6, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
                0, ...                                          % Initial phase shift
                true ...                                        % Polyphase implementation
                );
        case 'FBMC-QAM'
            Multicarrier = Modulation.FBMC(...
                3,...                                           % Number subcarriers.
                3,...                                           % Number FBMC symbols.
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmit real valued signal
                [PrototypeFilter '-QAM'],...                    % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
                3, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
                0, ...                                          % Initial phase shift
                true ...                                        % Polyphase implementation
                );
        case 'WOLA'
             Multicarrier = Modulation.WOLA(...
                floor(SamplingRate/SubcarrierSpacing/2),...     % Number subcarriers
                3,...                                           % Number WOLA Symbols
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmit real valued signal
                0, ...                                          % Cyclic prefix length (s) 
                0, ...                                          % Zero guard length (s)
                1/14/SubcarrierSpacing/2, ...                   % Length of the window overlapping (s) at the transmitter 
                1/14/SubcarrierSpacing/2 ...                    % Length of the window overlapping (s) at the receiver
                );
        case 'FOFDM'
             Multicarrier = Modulation.FOFDM(...
                floor(SamplingRate/SubcarrierSpacing/2),...     % Number subcarriers
                3,...                                           % Number FOFDM Symbols
                SubcarrierSpacing,...                           % Subcarrier spacing (Hz)
                SamplingRate,...                                % Sampling rate (Samples/s)
                0,...                                           % Intermediate frequency first subcarrier (Hz)
                false,...                                       % Transmit real valued signal
                0, ...                                          % Cyclic prefix length (s) 
                0, ...                                          % Zero guard length (s)
                1/14/SubcarrierSpacing/2, ...                   % Length of the transmit filter (s)
                1/14/SubcarrierSpacing/2, ...                   % Length of the receive filter (s)  
                1/14/SubcarrierSpacing ...                      % Length of the additional cyclic prefix (s). 
                );            
        otherwise
            error('Not supportet');
    end
    M_RealSubcarrierSpacing(i_subcarrier) = Multicarrier.PHY.SubcarrierSpacing;

    %% Channel Object
    ChannelModel = Channel.FastFading(...
        SamplingRate,...                                        % Sampling rate (Samples/s)
        PowerDelayProfile,...                                   % Power delay profile, either string or vector: 'Flat', 'AWGN', 'PedestrianA', 'PedestrianB', 'VehicularA', 'VehicularB', 'ExtendedPedestrianA', 'ExtendedPedestrianB', or 'TDL-A_xxns','TDL-B_xxns','TDL-C_xxns' (with xx the RMS delay spread in ns, e.g. 'TDL-A_30ns'), or [1 0 0.2] (Self-defined power delay profile which depends on the sampling rate)
        Multicarrier.Nr.SamplesTotal,...                        % Number of total samples
        Velocity_kmh/3.6*CarrierFrequency/2.998e8,...           % Maximum Doppler shift: Velocity_kmh/3.6*CarrierFrequency/2.998e8
        'Jakes',...                                             % Which Doppler model: 'Jakes', 'Uniform', 'Discrete-Jakes', 'Discrete-Uniform'. For "Discrete-", we assume a discrete Doppler spectrum to improve the simulation time. This only works accuratly if the number of samples and the velocity is sufficiently large
        200, ...                                                % Number multipath delays for WSS process
        1,...                                                   % Number of transmit antennas
        1,...                                                   % Number of receive antennas
        1 ...                                                   % Gives a warning if the predefined delay taps of the channel do not fit the sampling rate. This is usually not much of a problem if they are approximatly the same.
    );
    R_vecH = ChannelModel.GetCorrelationMatrix;
    TimeOffset = round(ChannelModel.GetMeanDelay/ChannelModel.PHY.dt);
%      TimeOffset = round(ChannelModel.GetRmsDelaySpread/ChannelModel.PHY.dt);
%      TimeOffset = 0;
    if  max(strcmp(Method,{'CP-OFDM','CP-OFDM-LTE','CP-OFDM-007LTE','CP-OFDM-010LTE'}))
        TimeOffset = 0;
    end
    
    switch Method
        case {'OFDM','CP-OFDM','FBMC-QAM','CP-OFDM-LTE','CP-OFDM-007LTE','CP-OFDM-010LTE','WOLA','FOFDM'}
            [PSignal_Theory,PInterference_Theory] = Multicarrier.GetSignalAndInterferencePowerQAM(...
                R_vecH,...
                eye(Multicarrier.Nr.Subcarriers*Multicarrier.Nr.MCSymbols),...
                TimeOffset,...
                ceil(Multicarrier.Nr.Subcarriers/2),...
                ceil(Multicarrier.Nr.MCSymbols/2));
            % Analytical solution for flat fading
            P_Signal_OFDM_Theory_Flat = hypergeom([1/2],[3/2 2],-(ChannelModel.PHY.MaximumDopplerShift/SubcarrierSpacing*pi).^2);
            P_Interference_OFDM_Theory_Flat = 1-P_Signal_OFDM_Theory_Flat;
            SIR_Theory_OFDM_Flat = 10*log10(P_Signal_OFDM_Theory_Flat/P_Interference_OFDM_Theory_Flat);
    case 'FBMC-OQAM'
            [PSignal_Theory,PInterference_Theory] = Multicarrier.GetSignalAndInterferencePowerOQAM(...
                R_vecH,...
                eye(Multicarrier.Nr.Subcarriers*Multicarrier.Nr.MCSymbols),...
                TimeOffset,...
                ceil(Multicarrier.Nr.Subcarriers/2),...
                ceil(Multicarrier.Nr.MCSymbols/2));
    end
    SIR_Theory(i_subcarrier) = 10*log10(PSignal_Theory./PInterference_Theory);
    disp(i_subcarrier/length(M_SubcarrierSpacing));
end

plot(M_RealSubcarrierSpacing/1e3,SIR_Theory); hold on;
xlabel('Subcarrier spacing (kHz)');
ylabel('Signal-to-Interference Ratio');
pause(0.1);

if SaveStuff
    FileName = ['./OptimalSubcarrierSpacing/Results/' Method '_' PrototypeFilter '_' PowerDelayProfile int2str(Velocity_kmh) '.mat'];
    save(FileName,...
          'SIR_Theory', ...
          'M_RealSubcarrierSpacing',...
          'Velocity_kmh', ...
          'CarrierFrequency', ...
          'M_SubcarrierSpacing',...
          'SamplingRate',...
          'Method',...
          'PrototypeFilter');  
end

end

