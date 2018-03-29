% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script calculates the Signal-to-Interference Ratio in case that two
% users share the same band. Both users utilize different subarrier
% spacings, 15kHz and 120kHz, to account for different channel conditions
% or because a low latency transmission is desired (120kHz). The SIR then
% depends on the guard band between those two users. We consider different 
% multicarrier modulation techniques:
%       1) OFDM with CP (worst spectral behaviour)
%       2) FBMC (best spectral properties, complex orthogonality is replaced by real orthogonality)
%       3) WOLA (windowed OFDM. The windowing is done at TX and RX)
%       4) FOFDM (filtered OFDM, sinc + Hann, filtering at TX and RX)
%       5) UFMC (filtered OFDM, subband-wise filtering, Dolph-Chebyshev window, cyclic prefix and zero padding are supported, filtering at TX and RX)
% See Figure 10 and 11 in the paper. For the case that no receive windowing
% and filtering is used, we have to set 
%   ReceiveFilterForWOLA_UFMC_FOFDM=false;

clear; close all;


M_NormalizedGuardBand                   = [0:0.025:0.5]; % [0:0.005:0.5];   % The normalized guard band F_G /(FL), wheras F_G represents the guard band between two use cases and FL the bandidth of each use case (assumed to be the same)
NormalizedGuardBandForWhichPSDisPlotted = 0.2;                              % For this value of the normalized guard band, the power spectral density is also calculated and plotted (presentation purpose only)

ReceiveFilterForWOLA_UFMC_FOFDM         = true;                             % Set to false to investigate effect if no receive filter is employed. Note that the TF-spacing could then be reduced (not done at the moment).
Bandwidth                               = '1.44MHz';                        % '1.44MHz' or '10.08MHz': load predefined values
PrototypeFilter                         = 'PHYDYAS';                        % PHYDYAS, Hermite, RRC
UFMC_ZeroPadding                        = false;                            % If true, apply zero padding instead of a cyclic prefix.
OverlappingFactorFBMC                   = 4;                                % Overlapping factor of FBMC. High overlapping improves the spectral behavior

if strcmp(Bandwidth,'1.44MHz')
    % Use Case 1: 
    L_15kHz                             = 12*8;                             % Number of subcarriers for use case 1
    K_15kHz                             = 9;                                % Number of time symbols for use case 1
    SubcarrierSpacing_15kHz             = 15e3; 

    % Use Case 2: 
    L_120kHz                            = L_15kHz/8;                        % Number of subcarriers for use case 2
    K_120kHz                            = K_15kHz*8;                        % Number of time symbols for use case 2
    SubcarrierSpacing_120kHz            = SubcarrierSpacing_15kHz*8;

    % Use Case 3 for FBMC. Note that 480kHz represents only the case for O=4
    L_480kHz                            = L_120kHz/OverlappingFactorFBMC;
    K_480kHz                            = K_120kHz*OverlappingFactorFBMC;
    SubcarrierSpacing_480kHz            = 120e3*OverlappingFactorFBMC;

    SamplingRate                        = SubcarrierSpacing_15kHz*L_15kHz*14;
    dt                                  = 1/SamplingRate;

    % Parameterset 1
    TF_FilteredAndWindowedOFDM_15kHz    = 1.09;                                    % Time frequency spacing for UFMC and FOFDM. Chosen so that the interference from filtering is small (the SIR is approximately 65dB).
    FilterLengthTXandRX_UFMC_15kHz      = 1/(14*SubcarrierSpacing_15kHz);          % 4.7탎, same as the CP in LTE. However, we need filtering at TX and RX => Total length of the filtering is two times larger
    FilterLengthTXandRX_FOFDM_15kHz     = 1.4*1/(14*SubcarrierSpacing_15kHz);      % Increased by 1.2 because FOFDM is more robust than UFMC to ISI (due to the full band filtering)

    % Parameterset 2
    TF_FilteredAndWindowedOFDM_120kHz   = 1.27;                                    % Time frequency spacing for UFMC and FOFDM. Chosen so that the interference from filtering is small (the SIR is approximately 65dB).
    FilterLengthTXandRX_UFMC_120kHz     = 3/(14*SubcarrierSpacing_120kHz);         % 3 times as high as the CP in LTE to get a smaller OOB. 
    FilterLengthTXandRX_FOFDM_120kHz    = 1.25*3/(14*SubcarrierSpacing_120kHz);    % Increased by 1.25 because FOFDM is more robust than UFMC to ISI (due to the filter)
elseif  strcmp(Bandwidth,'10.08MHz')
    % Use Case 1: 
    L_15kHz                             = 12*8*7;
    K_15kHz                             = 9;
    SubcarrierSpacing_15kHz             = 15e3;

    % Use Case 2: 
    L_120kHz                            = L_15kHz/8;
    K_120kHz                            = K_15kHz*8;
    SubcarrierSpacing_120kHz            = SubcarrierSpacing_15kHz*8;

    % Use Case 3 for FBMC. Note that 480kHz represents only the case for O=4
    L_480kHz                            = L_120kHz/OverlappingFactorFBMC;
    K_480kHz                            = K_120kHz*OverlappingFactorFBMC;
    SubcarrierSpacing_480kHz            = 120e3*OverlappingFactorFBMC;

    SamplingRate                        = SubcarrierSpacing_15kHz*L_15kHz*4;
    dt                                  = 1/SamplingRate;

    % Parameterset 1
    TF_FilteredAndWindowedOFDM_15kHz    = 1.09;                                     % Chosen so that the interference from the filtering is small enough (the SIR is approximately 65dB).
    FilterLengthTXandRX_UFMC_15kHz      = 1/(14*SubcarrierSpacing_15kHz);           % 4.7탎, same as the CP in LTE. But we need filtering at TX and RX => two times in total!
    FilterLengthTXandRX_FOFDM_15kHz     = 1.75*1/(14*SubcarrierSpacing_15kHz);      % Increased by 1.75 because FOFDM is more robust than UFMC to a larger filter length

    % Parameterset 2
    TF_FilteredAndWindowedOFDM_120kHz   = 1.09;                                     % Chosen so that the interference from the filtering is small enough (the SIR is approximately 65dB).
    FilterLengthTXandRX_UFMC_120kHz     = 1/(14*SubcarrierSpacing_120kHz);          % 3 times as high as the CP in LTE to get a smaller OOB. 
    FilterLengthTXandRX_FOFDM_120kHz    = 1.2*1/(14*SubcarrierSpacing_120kHz);      % Increased by 1.2 because FOFDM is more robust than UFMC to a larger filter length
else
    error('Not supported');
end
    
KL                                      = L_15kHz*K_15kHz;
CP_LengthFilter_15kHz                   = (TF_FilteredAndWindowedOFDM_15kHz-1)*1/(SubcarrierSpacing_15kHz);
CP_LengthFilter_120kHz                  = (TF_FilteredAndWindowedOFDM_120kHz-1)*1/(SubcarrierSpacing_120kHz);

ColorFBMC   = [0 0 1]*0.5;
ColorOFDM   = [1 0 0];
ColorWOLA   = [1 0 1];
ColorFOFDM  = [1 1 0]*0.7;
ColorUFMC   = [0 1 0]*0.5;

%% FBMC Object
FBMC_15kHz = Modulation.FBMC(...
    L_15kHz,...                                                              % Number subcarriers
    K_15kHz*2,...                                                            % Number FBMC symbols
    SubcarrierSpacing_15kHz,...                                              % Subcarrier spacing (Hz)
    SamplingRate,...                                                         % Sampling rate (Samples/s)
    0,...                                                                    % Intermediate frequency first subcarrier (Hz)
    false,...                                                                % Transmit real valued signal
    [PrototypeFilter '-OQAM'],...                                            % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    OverlappingFactorFBMC, ...                                               % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                                   % Initial phase shift
    true ...                                                                 % Polyphase implementation
    );
FBMC_120kHz = Modulation.FBMC(...
    L_120kHz,...                                                            % Number subcarriers
    K_120kHz*2,...                                                          % Number FBMC symbols
    SubcarrierSpacing_120kHz,...                                            % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    [PrototypeFilter '-OQAM'],...                                           % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    OverlappingFactorFBMC, ...                                              % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                                  % Initial phase shift
    true ...                                                                % Polyphase implementation
    );
FBMC_480kHz = Modulation.FBMC(...
    L_480kHz,...                                                            % Number subcarriers
    K_480kHz*2,...                                                          % Number FBMC symbols
    SubcarrierSpacing_480kHz,...                                            % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    [PrototypeFilter '-OQAM'],...                                           % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    OverlappingFactorFBMC, ...                                              % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                                  % Initial phase shift
    true ...                                                                % Polyphase implementation
    );

%% Conventional OFDM Object
OFDM_15kHz = Modulation.OFDM(...
    L_15kHz,...                                                             % Number subcarriers
    K_15kHz,...                                                             % Number OFDM Symbols
    SubcarrierSpacing_15kHz,...                                             % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    1/(14*SubcarrierSpacing_15kHz), ...                                     % Cyclic prefix length (s) 
    0 ...                                                                   % Zero guard length (s)
    );
OFDM_120kHz = Modulation.OFDM(...
    L_120kHz,...                                                            % Number subcarriers
    K_120kHz,...                                                            % Number OFDM Symbols
    SubcarrierSpacing_120kHz,...                                            % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    1/(14*SubcarrierSpacing_120kHz), ...                                    % Cyclic prefix length (s) 
    0 ...                                                                   % Zero guard length (s)
    );

%% Filtered OFDM Object
FOFDM_15kHz = Modulation.FOFDM(...
    L_15kHz,...                                                             % Number subcarriers
    K_15kHz,...                                                             % Number FOFDM Symbols
    SubcarrierSpacing_15kHz,...                                             % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    0, ...                                                                  % Cyclic prefix length (s) 
    0, ...                                                                  % Zero guard length (s)
    FilterLengthTXandRX_FOFDM_15kHz, ...                                    % Length of the transmit filter (s)
    ReceiveFilterForWOLA_UFMC_FOFDM*FilterLengthTXandRX_FOFDM_15kHz, ...    % Length of the receive filter (s) 
    CP_LengthFilter_15kHz...                                                % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
);
FOFDM_120kHz = Modulation.FOFDM(...
    L_120kHz,...                                                            % Number subcarriers
    K_120kHz,...                                                            % Number FOFDM Symbols
    SubcarrierSpacing_120kHz,...                                            % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    0, ...                                                                  % Cyclic prefix length (s) 
    0, ...                                                                  % Zero guard length (s)
    FilterLengthTXandRX_FOFDM_120kHz, ...                                   % Length of the transmit filter (s)    
    ReceiveFilterForWOLA_UFMC_FOFDM*FilterLengthTXandRX_FOFDM_120kHz, ...   % Length of the receive filter (s) 
    CP_LengthFilter_120kHz ...                                              % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
);

%% UFMC Oject
UFMC_15kHz = Modulation.UFMC(...
    L_15kHz,...                                                             % Number subcarriers
    K_15kHz,...                                                             % Number UFMC Symbols
    SubcarrierSpacing_15kHz,...                                             % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    0, ...                                                                  % Cyclic prefix length (s) 
    0, ...                                                                  % Zero guard length (s)
    FilterLengthTXandRX_UFMC_15kHz, ...                                     % Length of the transmit filter (s)
    ReceiveFilterForWOLA_UFMC_FOFDM*FilterLengthTXandRX_UFMC_15kHz, ...     % Length of the receive filter (s) 
    CP_LengthFilter_15kHz, ...                                              % Length of the additional cyclic prefix (or zero guard symbol if ZP is used) in seconds (s). Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
    UFMC_ZeroPadding ...                                                    % TRUE for Zero Padding (ZP) and FALSE for a conventional Cyclic Prefix (CP)
);
UFMC_120kHz = Modulation.UFMC(...
    L_120kHz,...                                                            % Number subcarriers
    K_120kHz,...                                                            % Number UFMC Symbols
    SubcarrierSpacing_120kHz,...                                            % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    0, ...                                                                  % Cyclic prefix length (s) 
    0, ...                                                                  % Zero guard length (s)
    FilterLengthTXandRX_UFMC_120kHz, ...                                    % Length of the transmit filter (s)
    ReceiveFilterForWOLA_UFMC_FOFDM*FilterLengthTXandRX_UFMC_120kHz, ...    % Length of the receive filter (s) 
    CP_LengthFilter_120kHz, ...                                             % Length of the additional cyclic prefix (or zero guard symbol if ZP is used) in seconds (s). Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
    UFMC_ZeroPadding ...                                                    % TRUE for Zero Padding (ZP) and FALSE for a conventional Cyclic Prefix (CP)
);


%% WOLA object
WOLA_15kHz = Modulation.WOLA(...
    L_15kHz,...                                                             % Number of subcarriers
    K_15kHz,...                                                             % Number WOLA Symbols
    SubcarrierSpacing_15kHz,...                                             % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    0, ...                                                                  % Cyclic prefix length (s) 
    0, ...                                                                  % Zero guard length (s)
    CP_LengthFilter_15kHz/2, ...                                            % Length of the window overlapping (s) at the transmitter 
    ReceiveFilterForWOLA_UFMC_FOFDM*CP_LengthFilter_15kHz/2 ...             % Length of the window overlapping (s) at the receiver
    );

WOLA_120kHz = Modulation.WOLA(...
    L_120kHz,...                                                            % Number subcarriers
    K_120kHz,...                                                            % Number WOLA Symbols
    SubcarrierSpacing_120kHz,...                                            % Subcarrier spacing (Hz)
    SamplingRate,...                                                        % Sampling rate (Samples/s)
    0,...                                                                   % Intermediate frequency first subcarrier (Hz)
    false,...                                                               % Transmit real valued signal
    0, ...                                                                  % Cyclic prefix length (s) 
    0, ...                                                                  % Zero guard length (s)
    CP_LengthFilter_120kHz/2, ...                                           % Length of the window overlapping (s) at the transmitter 
    ReceiveFilterForWOLA_UFMC_FOFDM*CP_LengthFilter_120kHz/2 ...            % Length of the window overlapping (s) at the receiver
    );


%% Plot Information
fprintf('==================================================================\n');
fprintf('=========================== 15kHz ================================\n');
fprintf('==================================================================\n');
fprintf('               |(complex)TF-Spacing| Bandwidth(LF)|   Time (KT)   \n');
fprintf('OFDM (with CP) |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', OFDM_15kHz.PHY.TimeSpacing*OFDM_15kHz.PHY.SubcarrierSpacing , OFDM_15kHz.PHY.SubcarrierSpacing*OFDM_15kHz.Nr.Subcarriers/1e6,OFDM_15kHz.PHY.TimeSpacing*OFDM_15kHz.Nr.MCSymbols/1e-6);
fprintf('FBMC           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', FBMC_15kHz.PHY.TimeSpacing*FBMC_15kHz.PHY.SubcarrierSpacing*2 , FBMC_15kHz.PHY.SubcarrierSpacing*FBMC_15kHz.Nr.Subcarriers/1e6,FBMC_15kHz.PHY.TimeSpacing*FBMC_15kHz.Nr.MCSymbols/1e-6);
fprintf('WOLA           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', WOLA_15kHz.PHY.TimeSpacing*WOLA_15kHz.PHY.SubcarrierSpacing , WOLA_15kHz.PHY.SubcarrierSpacing*WOLA_15kHz.Nr.Subcarriers/1e6,WOLA_15kHz.PHY.TimeSpacing*WOLA_15kHz.Nr.MCSymbols/1e-6);
fprintf('FOFDM          |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', FOFDM_15kHz.PHY.TimeSpacing*FOFDM_15kHz.PHY.SubcarrierSpacing , FOFDM_15kHz.PHY.SubcarrierSpacing*FOFDM_15kHz.Nr.Subcarriers/1e6,FOFDM_15kHz.PHY.TimeSpacing*FOFDM_15kHz.Nr.MCSymbols/1e-6);
fprintf('UFMC           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', UFMC_15kHz.PHY.TimeSpacing*UFMC_15kHz.PHY.SubcarrierSpacing , UFMC_15kHz.PHY.SubcarrierSpacing*UFMC_15kHz.Nr.Subcarriers/1e6,UFMC_15kHz.PHY.TimeSpacing*UFMC_15kHz.Nr.MCSymbols/1e-6);

fprintf('==================================================================\n');
fprintf('=========================== 120kHz ===============================\n');
fprintf('==================================================================\n');
fprintf('               |(complex)TF-Spacing| Bandwidth(LF)|   Time (KT) \n');
fprintf('OFDM (with CP) |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', OFDM_120kHz.PHY.TimeSpacing*OFDM_120kHz.PHY.SubcarrierSpacing , OFDM_120kHz.PHY.SubcarrierSpacing*OFDM_120kHz.Nr.Subcarriers/1e6,OFDM_120kHz.PHY.TimeSpacing*OFDM_120kHz.Nr.MCSymbols/1e-6);
fprintf('FBMC           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', FBMC_120kHz.PHY.TimeSpacing*FBMC_120kHz.PHY.SubcarrierSpacing*2 , FBMC_120kHz.PHY.SubcarrierSpacing*FBMC_120kHz.Nr.Subcarriers/1e6,FBMC_120kHz.PHY.TimeSpacing*FBMC_120kHz.Nr.MCSymbols/1e-6);
fprintf('WOLA           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', WOLA_120kHz.PHY.TimeSpacing*WOLA_120kHz.PHY.SubcarrierSpacing , WOLA_120kHz.PHY.SubcarrierSpacing*WOLA_120kHz.Nr.Subcarriers/1e6,WOLA_120kHz.PHY.TimeSpacing*WOLA_120kHz.Nr.MCSymbols/1e-6);
fprintf('FOFDM          |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', FOFDM_120kHz.PHY.TimeSpacing*FOFDM_120kHz.PHY.SubcarrierSpacing , FOFDM_120kHz.PHY.SubcarrierSpacing*FOFDM_120kHz.Nr.Subcarriers/1e6,FOFDM_120kHz.PHY.TimeSpacing*FOFDM_120kHz.Nr.MCSymbols/1e-6);
fprintf('UFMC           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', UFMC_120kHz.PHY.TimeSpacing*UFMC_120kHz.PHY.SubcarrierSpacing , UFMC_120kHz.PHY.SubcarrierSpacing*UFMC_120kHz.Nr.Subcarriers/1e6,UFMC_120kHz.PHY.TimeSpacing*UFMC_120kHz.Nr.MCSymbols/1e-6);

fprintf('==================================================================\n');
fprintf('=========================== 480kHz ===============================\n');
fprintf('==================================================================\n');
fprintf('               |(complex)TF-Spacing| Bandwidth(LF)|   Time (KT) \n');
fprintf('FBMC           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', FBMC_480kHz.PHY.TimeSpacing*FBMC_480kHz.PHY.SubcarrierSpacing*2 , FBMC_480kHz.PHY.SubcarrierSpacing*FBMC_480kHz.Nr.Subcarriers/1e6,FBMC_480kHz.PHY.TimeSpacing*FBMC_480kHz.Nr.MCSymbols/1e-6);


%% Generate Transmit  Matrices
G_TX_OFDM_15kHz_0Hz_DiffLength  = sparse(OFDM_15kHz.GetTXMatrix);
G_TX_FBMC_15kHz_0Hz_DiffLength  = sparse(FBMC_15kHz.GetTXMatrix);
G_TX_WOLA_15kHz_0Hz_DiffLength  = sparse(WOLA_15kHz.GetTXMatrix);
G_TX_FOFDM_15kHz_0Hz_DiffLength = sparse(FOFDM_15kHz.GetTXMatrix);
G_TX_UFMC_15kHz_0Hz_DiffLength  = sparse(UFMC_15kHz.GetTXMatrix);

G_TX_OFDM_120kHz_0Hz_DiffLength = sparse(OFDM_120kHz.GetTXMatrix);
G_TX_FBMC_120kHz_0Hz_DiffLength = sparse(FBMC_120kHz.GetTXMatrix);
G_TX_WOLA_120kHz_0Hz_DiffLength = sparse(WOLA_120kHz.GetTXMatrix);
G_TX_FOFDM_120kHz_0Hz_DiffLength= sparse(FOFDM_120kHz.GetTXMatrix);
G_TX_UFMC_120kHz_0Hz_DiffLength = sparse(UFMC_120kHz.GetTXMatrix);

G_TX_FBMC_480kHz_0Hz_DiffLength = sparse(FBMC_480kHz.GetTXMatrix);


disp('Generate receive matrices ... this may take a while...');
G_RX_OFDM_15kHz_0Hz_DiffLength  = sparse(OFDM_15kHz.GetRXMatrix);
G_RX_FBMC_15kHz_0Hz_DiffLength  = sparse(FBMC_15kHz.GetRXMatrix);
G_RX_WOLA_15kHz_0Hz_DiffLength  = sparse(WOLA_15kHz.GetRXMatrix);
G_RX_FOFDM_15kHz_0Hz_DiffLength = sparse(FOFDM_15kHz.GetRXMatrix);
G_RX_UFMC_15kHz_0Hz_DiffLength  = sparse(UFMC_15kHz.GetRXMatrix);

G_RX_OFDM_120kHz_0Hz_DiffLength = sparse(OFDM_120kHz.GetRXMatrix);
G_RX_FBMC_120kHz_0Hz_DiffLength = sparse(FBMC_120kHz.GetRXMatrix);
G_RX_WOLA_120kHz_0Hz_DiffLength = sparse(WOLA_120kHz.GetRXMatrix);
G_RX_FOFDM_120kHz_0Hz_DiffLength= sparse(FOFDM_120kHz.GetRXMatrix);
G_RX_UFMC_120kHz_0Hz_DiffLength = sparse(UFMC_120kHz.GetRXMatrix);

G_RX_FBMC_480kHz_0Hz_DiffLength = sparse(FBMC_480kHz.GetRXMatrix);


% Add zeros so that 15kHz and 120kHz fit together in time
if size(G_TX_OFDM_15kHz_0Hz_DiffLength,1)>size(G_TX_OFDM_120kHz_0Hz_DiffLength,1)
    Ndiff                   = size(G_TX_OFDM_15kHz_0Hz_DiffLength,1)-size(G_TX_OFDM_120kHz_0Hz_DiffLength,1);
    
    G_TX_OFDM_15kHz_0Hz     = G_TX_OFDM_15kHz_0Hz_DiffLength;
    G_TX_OFDM_120kHz_0Hz    = [zeros(ceil(Ndiff/2),KL);G_TX_OFDM_120kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL)];
    
    G_RX_OFDM_15kHz_0Hz     = G_RX_OFDM_15kHz_0Hz_DiffLength;
    G_RX_OFDM_120kHz_0Hz    = [zeros(ceil(Ndiff/2),KL)',G_RX_OFDM_120kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL)'];    
else
    Ndiff                   = size(G_TX_OFDM_120kHz_0Hz_DiffLength,1)-size(G_TX_OFDM_15kHz_0Hz_DiffLength,1);
    G_TX_OFDM_120kHz_0Hz    = G_TX_OFDM_120kHz_0Hz_DiffLength; 
    G_TX_OFDM_15kHz_0Hz     = [zeros(ceil(Ndiff/2),KL);G_TX_OFDM_15kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL)];
    
    G_RX_OFDM_120kHz_0Hz    = G_RX_OFDM_120kHz_0Hz_DiffLength; 
    G_RX_OFDM_15kHz_0Hz     = [zeros(ceil(Ndiff/2),KL)',G_RX_OFDM_15kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL)'];   
end
if size(G_TX_FBMC_15kHz_0Hz_DiffLength,1)>size(G_TX_FBMC_120kHz_0Hz_DiffLength,1)
    Ndiff                   = size(G_TX_FBMC_15kHz_0Hz_DiffLength,1)-size(G_TX_FBMC_120kHz_0Hz_DiffLength,1);
    
    G_TX_FBMC_15kHz_0Hz     = G_TX_FBMC_15kHz_0Hz_DiffLength;
    G_TX_FBMC_120kHz_0Hz    = [zeros(ceil(Ndiff/2),KL*2);G_TX_FBMC_120kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL*2)];
    
    G_RX_FBMC_15kHz_0Hz     = G_RX_FBMC_15kHz_0Hz_DiffLength;
    G_RX_FBMC_120kHz_0Hz    = [zeros(ceil(Ndiff/2),KL*2)',G_RX_FBMC_120kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL*2)'];    
else
    Ndiff                   = size(G_TX_FBMC_120kHz_0Hz_DiffLength,1)-size(G_TX_FBMC_15kHz_0Hz_DiffLength,1);
    G_TX_FBMC_120kHz_0Hz    = G_TX_FBMC_120kHz_0Hz_DiffLength; 
    G_TX_FBMC_15kHz_0Hz     = [zeros(ceil(Ndiff/2),KL*2);G_TX_FBMC_15kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL*2)];
    
    G_RX_FBMC_120kHz_0Hz    = G_RX_FBMC_120kHz_0Hz_DiffLength; 
    G_RX_FBMC_15kHz_0Hz     = [zeros(ceil(Ndiff/2),KL*2)',G_RX_FBMC_15kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL*2)'];   
end
if size(G_TX_WOLA_15kHz_0Hz_DiffLength,1)>size(G_TX_WOLA_120kHz_0Hz_DiffLength,1)
    Ndiff                   = size(G_TX_WOLA_15kHz_0Hz_DiffLength,1)-size(G_TX_WOLA_120kHz_0Hz_DiffLength,1);
    
    G_TX_WOLA_15kHz_0Hz     = G_TX_WOLA_15kHz_0Hz_DiffLength;
    G_TX_WOLA_120kHz_0Hz    = [zeros(ceil(Ndiff/2),KL);G_TX_WOLA_120kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL)];
    
    G_RX_WOLA_15kHz_0Hz     = G_RX_WOLA_15kHz_0Hz_DiffLength;
    G_RX_WOLA_120kHz_0Hz    = [zeros(ceil(Ndiff/2),KL)',G_RX_WOLA_120kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL)'];    
else
    Ndiff                   = size(G_TX_WOLA_120kHz_0Hz_DiffLength,1)-size(G_TX_WOLA_15kHz_0Hz_DiffLength,1);
    G_TX_WOLA_120kHz_0Hz    = G_TX_WOLA_120kHz_0Hz_DiffLength; 
    G_TX_WOLA_15kHz_0Hz     = [zeros(ceil(Ndiff/2),KL);G_TX_WOLA_15kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL)];
    
    G_RX_WOLA_120kHz_0Hz    = G_RX_WOLA_120kHz_0Hz_DiffLength; 
    G_RX_WOLA_15kHz_0Hz     = [zeros(ceil(Ndiff/2),KL)',G_RX_WOLA_15kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL)'];   
end
if size(G_TX_FOFDM_15kHz_0Hz_DiffLength,1)>size(G_TX_FOFDM_120kHz_0Hz_DiffLength,1)
    Ndiff                   = size(G_TX_FOFDM_15kHz_0Hz_DiffLength,1)-size(G_TX_FOFDM_120kHz_0Hz_DiffLength,1);
    
    G_TX_FOFDM_15kHz_0Hz    = G_TX_FOFDM_15kHz_0Hz_DiffLength;
    G_TX_FOFDM_120kHz_0Hz   = [zeros(ceil(Ndiff/2),KL);G_TX_FOFDM_120kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL)];
    
    G_RX_FOFDM_15kHz_0Hz    = G_RX_FOFDM_15kHz_0Hz_DiffLength;
    G_RX_FOFDM_120kHz_0Hz   = [zeros(ceil(Ndiff/2),KL)',G_RX_FOFDM_120kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL)'];    
else
    Ndiff                   = size(G_TX_FOFDM_120kHz_0Hz_DiffLength,1)-size(G_TX_FOFDM_15kHz_0Hz_DiffLength,1);
    G_TX_FOFDM_120kHz_0Hz   = G_TX_FOFDM_120kHz_0Hz_DiffLength; 
    G_TX_FOFDM_15kHz_0Hz    = [zeros(ceil(Ndiff/2),KL);G_TX_FOFDM_15kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL)];
    
    G_RX_FOFDM_120kHz_0Hz   = G_RX_FOFDM_120kHz_0Hz_DiffLength; 
    G_RX_FOFDM_15kHz_0Hz    = [zeros(ceil(Ndiff/2),KL)',G_RX_FOFDM_15kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL)'];   
end
if size(G_TX_UFMC_15kHz_0Hz_DiffLength,1)>size(G_TX_UFMC_120kHz_0Hz_DiffLength,1)
    Ndiff                   = size(G_TX_UFMC_15kHz_0Hz_DiffLength,1)-size(G_TX_UFMC_120kHz_0Hz_DiffLength,1);
    
    G_TX_UFMC_15kHz_0Hz     = G_TX_UFMC_15kHz_0Hz_DiffLength;
    G_TX_UFMC_120kHz_0Hz    = [zeros(ceil(Ndiff/2),KL);G_TX_UFMC_120kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL)];
    
    G_RX_UFMC_15kHz_0Hz     = G_RX_UFMC_15kHz_0Hz_DiffLength;
    G_RX_UFMC_120kHz_0Hz    = [zeros(ceil(Ndiff/2),KL)',G_RX_UFMC_120kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL)'];    
else
    Ndiff                   = size(G_TX_UFMC_120kHz_0Hz_DiffLength,1)-size(G_TX_UFMC_15kHz_0Hz_DiffLength,1);
    G_TX_UFMC_120kHz_0Hz    = G_TX_UFMC_120kHz_0Hz_DiffLength; 
    G_TX_UFMC_15kHz_0Hz     = [zeros(ceil(Ndiff/2),KL);G_TX_UFMC_15kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL)];
    
    G_RX_UFMC_120kHz_0Hz    = G_RX_UFMC_120kHz_0Hz_DiffLength; 
    G_RX_UFMC_15kHz_0Hz     = [zeros(ceil(Ndiff/2),KL)',G_RX_UFMC_15kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL)'];   
end

Ndiff                       = size(G_TX_FBMC_15kHz_0Hz,1)-size(G_TX_FBMC_480kHz_0Hz_DiffLength,1);
G_TX_FBMC_480kHz_0Hz        = [zeros(ceil(Ndiff/2),KL*2);G_TX_FBMC_480kHz_0Hz_DiffLength;zeros(floor(Ndiff/2),KL*2)];
G_RX_FBMC_480kHz_0Hz        = [zeros(ceil(Ndiff/2),KL*2)',G_RX_FBMC_480kHz_0Hz_DiffLength,zeros(floor(Ndiff/2),KL*2)']; 


%% The first subcarrier is at 0Hz. However, we have to perform a frequency shift by F/2 so that the whole bandwidth FL starts at 0Hz and not just the center of the first subcarrier
t_OFDM  =  (0:size(G_TX_OFDM_15kHz_0Hz,1)-1)*dt;
t_FBMC  =  (0:size(G_TX_FBMC_15kHz_0Hz,1)-1)*dt;
t_WOLA  =  (0:size(G_TX_WOLA_15kHz_0Hz,1)-1)*dt;
t_FOFDM =  (0:size(G_TX_FOFDM_15kHz_0Hz,1)-1)*dt;
t_UFMC  =  (0:size(G_TX_UFMC_15kHz_0Hz,1)-1)*dt;

G_TX_OFDM_15kHz     = bsxfun(@times,G_TX_OFDM_15kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_15kHz/2*t_OFDM.'));
G_RX_OFDM_15kHz     = bsxfun(@times,G_RX_OFDM_15kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_15kHz/2*t_OFDM));
G_TX_OFDM_120kHz    = bsxfun(@times,G_TX_OFDM_120kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_120kHz/2*t_OFDM.'));
G_RX_OFDM_120kHz    = bsxfun(@times,G_RX_OFDM_120kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_120kHz/2*t_OFDM));

G_TX_FBMC_15kHz     = bsxfun(@times,G_TX_FBMC_15kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_15kHz/2*t_FBMC.'));
G_RX_FBMC_15kHz     = bsxfun(@times,G_RX_FBMC_15kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_15kHz/2*t_FBMC));
G_TX_FBMC_120kHz    = bsxfun(@times,G_TX_FBMC_120kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_120kHz/2*t_FBMC.'));
G_RX_FBMC_120kHz    = bsxfun(@times,G_RX_FBMC_120kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_120kHz/2*t_FBMC));
G_TX_FBMC_480kHz    = bsxfun(@times,G_TX_FBMC_480kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_480kHz/2*t_FBMC.'));
G_RX_FBMC_480kHz    = bsxfun(@times,G_RX_FBMC_480kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_480kHz/2*t_FBMC));

G_TX_WOLA_15kHz     = bsxfun(@times,G_TX_WOLA_15kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_15kHz/2*t_WOLA.'));
G_RX_WOLA_15kHz     = bsxfun(@times,G_RX_WOLA_15kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_15kHz/2*t_WOLA));
G_TX_WOLA_120kHz    = bsxfun(@times,G_TX_WOLA_120kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_120kHz/2*t_WOLA.'));
G_RX_WOLA_120kHz    = bsxfun(@times,G_RX_WOLA_120kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_120kHz/2*t_WOLA));

G_TX_FOFDM_15kHz    = bsxfun(@times,G_TX_FOFDM_15kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_15kHz/2*t_FOFDM.'));
G_RX_FOFDM_15kHz    = bsxfun(@times,G_RX_FOFDM_15kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_15kHz/2*t_FOFDM));
G_TX_FOFDM_120kHz   = bsxfun(@times,G_TX_FOFDM_120kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_120kHz/2*t_FOFDM.'));
G_RX_FOFDM_120kHz   = bsxfun(@times,G_RX_FOFDM_120kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_120kHz/2*t_FOFDM));

G_TX_UFMC_15kHz     = bsxfun(@times,G_TX_UFMC_15kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_15kHz/2*t_UFMC.'));
G_RX_UFMC_15kHz     = bsxfun(@times,G_RX_UFMC_15kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_15kHz/2*t_UFMC));
G_TX_UFMC_120kHz    = bsxfun(@times,G_TX_UFMC_120kHz_0Hz,exp(1j*2*pi*SubcarrierSpacing_120kHz/2*t_UFMC.'));
G_RX_UFMC_120kHz    = bsxfun(@times,G_RX_UFMC_120kHz_0Hz,exp(-1j*2*pi*SubcarrierSpacing_120kHz/2*t_UFMC));



%% Calculate the in-band transmission matrix
D_OFDM_15kHz    = G_RX_OFDM_15kHz*G_TX_OFDM_15kHz;
D_OFDM_120kHz   = G_RX_OFDM_120kHz*G_TX_OFDM_120kHz;

D_FBMC_15kHz    = G_RX_FBMC_15kHz*G_TX_FBMC_15kHz;
D_FBMC_120kHz   = G_RX_FBMC_120kHz*G_TX_FBMC_120kHz;
D_FBMC_480kHz   = G_RX_FBMC_480kHz*G_TX_FBMC_480kHz;

D_WOLA_15kHz    = G_RX_WOLA_15kHz*G_TX_WOLA_15kHz;
D_WOLA_120kHz   = G_RX_WOLA_120kHz*G_TX_WOLA_120kHz;

D_FOFDM_15kHz   = G_RX_FOFDM_15kHz*G_TX_FOFDM_15kHz;
D_FOFDM_120kHz  = G_RX_FOFDM_120kHz*G_TX_FOFDM_120kHz;

D_UFMC_15kHz    = G_RX_UFMC_15kHz*G_TX_UFMC_15kHz;
D_UFMC_120kHz   = G_RX_UFMC_120kHz*G_TX_UFMC_120kHz;

%% Calculate the in-band SIR 
SignalPowerInBand_OFDM_15kHz        = full(reshape(abs(diag(D_OFDM_15kHz)).^2,L_15kHz,K_15kHz));
InterferencePowerInBand_OFDM_15kHz  = reshape(sum(abs(D_OFDM_15kHz).^2,2),L_15kHz,K_15kHz) - SignalPowerInBand_OFDM_15kHz;
SignalPowerInBand_OFDM_120kHz       = full(reshape(abs(diag(D_OFDM_120kHz)).^2,L_120kHz,K_120kHz));
InterferencePowerInBand_OFDM_120kHz = reshape(sum(abs(D_OFDM_120kHz).^2,2),L_120kHz,K_120kHz) - SignalPowerInBand_OFDM_120kHz;

SignalPowerInBand_FBMC_15kHz        = full(reshape(abs(diag(real(D_FBMC_15kHz))).^2,L_15kHz,2*K_15kHz));
InterferencePowerInBand_FBMC_15kHz  = reshape(sum(abs(real(D_FBMC_15kHz)).^2,2),L_15kHz,2*K_15kHz) - SignalPowerInBand_FBMC_15kHz;
SignalPowerInBand_FBMC_120kHz       = full(reshape(abs(diag(real(D_FBMC_120kHz))).^2,L_120kHz,2*K_120kHz));
InterferencePowerInBand_FBMC_120kHz = reshape(sum(abs(real(D_FBMC_120kHz)).^2,2),L_120kHz,2*K_120kHz) - SignalPowerInBand_FBMC_120kHz;
SignalPowerInBand_FBMC_480kHz       = full(reshape(abs(diag(real(D_FBMC_480kHz))).^2,L_480kHz,2*K_480kHz));
InterferencePowerInBand_FBMC_480kHz = reshape(sum(abs(real(D_FBMC_480kHz)).^2,2),L_480kHz,2*K_480kHz) - SignalPowerInBand_FBMC_480kHz;


SignalPowerInBand_WOLA_15kHz        = full(reshape(abs(diag(D_WOLA_15kHz)).^2,L_15kHz,K_15kHz));
InterferencePowerInBand_WOLA_15kHz  = reshape(sum(abs(D_WOLA_15kHz).^2,2),L_15kHz,K_15kHz) - SignalPowerInBand_WOLA_15kHz;
SignalPowerInBand_WOLA_120kHz       = full(reshape(abs(diag(D_WOLA_120kHz)).^2,L_120kHz,K_120kHz));
InterferencePowerInBand_WOLA_120kHz = reshape(sum(abs(D_WOLA_120kHz).^2,2),L_120kHz,K_120kHz) - SignalPowerInBand_WOLA_120kHz;

SignalPowerInBand_FOFDM_15kHz       = full(reshape(abs(diag(D_FOFDM_15kHz)).^2,L_15kHz,K_15kHz));
InterferencePowerInBand_FOFDM_15kHz = reshape(sum(abs(D_FOFDM_15kHz).^2,2),L_15kHz,K_15kHz) - SignalPowerInBand_FOFDM_15kHz;
SignalPowerInBand_FOFDM_120kHz      = full(reshape(abs(diag(D_FOFDM_120kHz)).^2,L_120kHz,K_120kHz));
InterferencePowerInBand_FOFDM_120kHz= reshape(sum(abs(D_FOFDM_120kHz).^2,2),L_120kHz,K_120kHz) - SignalPowerInBand_FOFDM_120kHz;

SignalPowerInBand_UFMC_15kHz        = full(reshape(abs(diag(D_UFMC_15kHz)).^2,L_15kHz,K_15kHz));
InterferencePowerInBand_UFMC_15kHz  = reshape(sum(abs(D_UFMC_15kHz).^2,2),L_15kHz,K_15kHz) - SignalPowerInBand_UFMC_15kHz;
SignalPowerInBand_UFMC_120kHz       = full(reshape(abs(diag(D_UFMC_120kHz)).^2,L_120kHz,K_120kHz));
InterferencePowerInBand_UFMC_120kHz = reshape(sum(abs(D_UFMC_120kHz).^2,2),L_120kHz,K_120kHz) - SignalPowerInBand_UFMC_120kHz;


fprintf('=============== In-Band SIR  ======================\n');
fprintf('         OFDM  |   FBMC  |  WOLA  | FOFDM  |  UFMC  | \n');
fprintf('15kHz:  %4.2fdB | %4.2fdB | %4.2fdB | %4.2fdB | %4.2f dB | \n', 10*log10(sum(SignalPowerInBand_OFDM_15kHz(:))./sum(InterferencePowerInBand_OFDM_15kHz(:))), 10*log10(sum(SignalPowerInBand_FBMC_15kHz(:))./sum(InterferencePowerInBand_FBMC_15kHz(:))), 10*log10(sum(SignalPowerInBand_WOLA_15kHz(:))./sum(InterferencePowerInBand_WOLA_15kHz(:))),10*log10(sum(SignalPowerInBand_FOFDM_15kHz(:))./sum(InterferencePowerInBand_FOFDM_15kHz(:))),10*log10(sum(SignalPowerInBand_UFMC_15kHz(:))./sum(InterferencePowerInBand_UFMC_15kHz(:))));
fprintf('120kHz: %4.2fdB | %4.2fdB | %4.2fdB | %4.2fdB | %4.2f dB | \n', 10*log10(sum(SignalPowerInBand_OFDM_120kHz(:))./sum(InterferencePowerInBand_OFDM_120kHz(:))), 10*log10(sum(SignalPowerInBand_FBMC_120kHz(:))./sum(InterferencePowerInBand_FBMC_120kHz(:))), 10*log10(sum(SignalPowerInBand_WOLA_120kHz(:))./sum(InterferencePowerInBand_WOLA_120kHz(:))),10*log10(sum(SignalPowerInBand_FOFDM_120kHz(:))./sum(InterferencePowerInBand_FOFDM_120kHz(:))),10*log10(sum(SignalPowerInBand_UFMC_120kHz(:))./sum(InterferencePowerInBand_UFMC_120kHz(:))));
fprintf('480kHz:        | %4.2fdB |         |         |          | \n', 10*log10(sum(SignalPowerInBand_FBMC_480kHz(:))./sum(InterferencePowerInBand_FBMC_480kHz(:))));


for i_GuardBand = 1:length(M_NormalizedGuardBand)

%% Calculate the Interference 
FrequencyShiftNextBand              = L_15kHz*SubcarrierSpacing_15kHz;
GuardBand                           = FrequencyShiftNextBand*M_NormalizedGuardBand(i_GuardBand); 


G_TX_OFDM_120kHz_FrequencyShifted   = bsxfun(@times,G_TX_OFDM_120kHz,exp(1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_OFDM.'));
G_RX_OFDM_120kHz_FrequencyShifted   = bsxfun(@times,G_RX_OFDM_120kHz,exp(-1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_OFDM));

G_TX_FBMC_120kHz_FrequencyShifted   = bsxfun(@times,G_TX_FBMC_120kHz,exp(1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_FBMC.'));
G_RX_FBMC_120kHz_FrequencyShifted   = bsxfun(@times,G_RX_FBMC_120kHz,exp(-1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_FBMC));
G_TX_FBMC_480kHz_FrequencyShifted   = bsxfun(@times,G_TX_FBMC_480kHz,exp(1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_FBMC.'));
G_RX_FBMC_480kHz_FrequencyShifted   = bsxfun(@times,G_RX_FBMC_480kHz,exp(-1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_FBMC));

G_TX_WOLA_120kHz_FrequencyShifted   = bsxfun(@times,G_TX_WOLA_120kHz,exp(1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_WOLA.'));
G_RX_WOLA_120kHz_FrequencyShifted   = bsxfun(@times,G_RX_WOLA_120kHz,exp(-1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_WOLA));

G_TX_FOFDM_120kHz_FrequencyShifted  = bsxfun(@times,G_TX_FOFDM_120kHz,exp(1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_FOFDM.'));
G_RX_FOFDM_120kHz_FrequencyShifted  = bsxfun(@times,G_RX_FOFDM_120kHz,exp(-1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_FOFDM));

G_TX_UFMC_120kHz_FrequencyShifted   = bsxfun(@times,G_TX_UFMC_120kHz,exp(1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_UFMC.'));
G_RX_UFMC_120kHz_FrequencyShifted   = bsxfun(@times,G_RX_UFMC_120kHz,exp(-1j*2*pi*(FrequencyShiftNextBand+GuardBand)*t_UFMC));


% Calculate Interference
InterferencePower_120kHzTo15kHz_OFDM    = reshape(sum(abs(G_RX_OFDM_15kHz*G_TX_OFDM_120kHz_FrequencyShifted).^2,2),L_15kHz,K_15kHz);
InterferencePower_120kHzTo15kHz_FBMC    = reshape(sum(abs(real(G_RX_FBMC_15kHz*G_TX_FBMC_120kHz_FrequencyShifted)).^2,2),L_15kHz,2*K_15kHz);
InterferencePower_120kHzTo15kHz_WOLA    = reshape(sum(abs((G_RX_WOLA_15kHz*G_TX_WOLA_120kHz_FrequencyShifted)).^2,2),L_15kHz,K_15kHz);
InterferencePower_120kHzTo15kHz_FOFDM   = reshape(sum(abs((G_RX_FOFDM_15kHz*G_TX_FOFDM_120kHz_FrequencyShifted)).^2,2),L_15kHz,K_15kHz);
InterferencePower_120kHzTo15kHz_UFMC    = reshape(sum(abs((G_RX_UFMC_15kHz*G_TX_UFMC_120kHz_FrequencyShifted)).^2,2),L_15kHz,K_15kHz);


InterferencePower_15kHzTo120kHz_OFDM    = reshape(sum(abs(G_RX_OFDM_120kHz_FrequencyShifted*G_TX_OFDM_15kHz).^2,2),L_120kHz,K_120kHz);
InterferencePower_15kHzTo120kHz_FBMC    = reshape(sum(abs(real(G_RX_FBMC_120kHz_FrequencyShifted*G_TX_FBMC_15kHz)).^2,2),L_120kHz,2*K_120kHz);
InterferencePower_15kHzTo120kHz_WOLA    = reshape(sum(abs((G_RX_WOLA_120kHz_FrequencyShifted*G_TX_WOLA_15kHz)).^2,2),L_120kHz,K_120kHz);
InterferencePower_15kHzTo120kHz_FOFDM   = reshape(sum(abs((G_RX_FOFDM_120kHz_FrequencyShifted*G_TX_FOFDM_15kHz)).^2,2),L_120kHz,K_120kHz);
InterferencePower_15kHzTo120kHz_UFMC    = reshape(sum(abs((G_RX_UFMC_120kHz_FrequencyShifted*G_TX_UFMC_15kHz)).^2,2),L_120kHz,K_120kHz);

InterferencePower_480kHzTo15kHz_FBMC    = reshape(sum(abs(real(G_RX_FBMC_15kHz*G_TX_FBMC_480kHz_FrequencyShifted)).^2,2),L_15kHz,2*K_15kHz);
InterferencePower_15kHzTo480kHz_FBMC    = reshape(sum(abs(real(G_RX_FBMC_480kHz_FrequencyShifted*G_TX_FBMC_15kHz)).^2,2),L_480kHz,2*K_480kHz);


% SIR_dB  for each Position
SIR_dB_AllPositions_OFDM_15kHz(:,:,i_GuardBand)     = 10*log10(SignalPowerInBand_OFDM_15kHz./(InterferencePower_120kHzTo15kHz_OFDM+InterferencePowerInBand_OFDM_15kHz));
SIR_dB_AllPositions_FBMC_15kHz(:,:,i_GuardBand)     = 10*log10(SignalPowerInBand_FBMC_15kHz./(InterferencePower_120kHzTo15kHz_FBMC+InterferencePowerInBand_FBMC_15kHz));
SIR_dB_AllPositions_WOLA_15kHz(:,:,i_GuardBand)     = 10*log10(SignalPowerInBand_WOLA_15kHz./(InterferencePower_120kHzTo15kHz_WOLA+InterferencePowerInBand_WOLA_15kHz));
SIR_dB_AllPositions_FOFDM_15kHz(:,:,i_GuardBand)    = 10*log10(SignalPowerInBand_FOFDM_15kHz./(InterferencePower_120kHzTo15kHz_FOFDM+InterferencePowerInBand_FOFDM_15kHz));
SIR_dB_AllPositions_UFMC_15kHz(:,:,i_GuardBand)     = 10*log10(SignalPowerInBand_UFMC_15kHz./(InterferencePower_120kHzTo15kHz_UFMC+InterferencePowerInBand_UFMC_15kHz));

SIR_dB_AllPositions_OFDM_120kHz(:,:,i_GuardBand)    = 10*log10(SignalPowerInBand_OFDM_120kHz./(InterferencePower_15kHzTo120kHz_OFDM+InterferencePowerInBand_OFDM_120kHz));
SIR_dB_AllPositions_FBMC_120kHz(:,:,i_GuardBand)    = 10*log10(SignalPowerInBand_FBMC_120kHz./(InterferencePower_15kHzTo120kHz_FBMC+InterferencePowerInBand_FBMC_120kHz));
SIR_dB_AllPositions_WOLA_120kHz(:,:,i_GuardBand)    = 10*log10(SignalPowerInBand_WOLA_120kHz./(InterferencePower_15kHzTo120kHz_WOLA+InterferencePowerInBand_WOLA_120kHz));
SIR_dB_AllPositions_FOFDM_120kHz(:,:,i_GuardBand)   = 10*log10(SignalPowerInBand_FOFDM_120kHz./(InterferencePower_15kHzTo120kHz_FOFDM+InterferencePowerInBand_FOFDM_120kHz));
SIR_dB_AllPositions_UFMC_120kHz(:,:,i_GuardBand)    = 10*log10(SignalPowerInBand_UFMC_120kHz./(InterferencePower_15kHzTo120kHz_UFMC+InterferencePowerInBand_UFMC_120kHz));

SIR_dB_AllPositions_FBMC_15kHzFrom480kHz(:,:,i_GuardBand) = 10*log10(SignalPowerInBand_FBMC_15kHz./(InterferencePower_480kHzTo15kHz_FBMC+InterferencePowerInBand_FBMC_15kHz));
SIR_dB_AllPositions_FBMC_480kHzFrom15kHz(:,:,i_GuardBand) = 10*log10(SignalPowerInBand_FBMC_480kHz./(InterferencePower_15kHzTo480kHz_FBMC+InterferencePowerInBand_FBMC_480kHz));

% Signal to Interference Ratio over Frequency (time is averaged out and we remove the first 2 and last 2 symbols as a precaution)
SIR_dB_Frequency_OFDM_15kHz(:,i_GuardBand)      = 10*log10(sum(SignalPowerInBand_OFDM_15kHz(:,3:end-2),2)./(sum(InterferencePower_120kHzTo15kHz_OFDM(:,3:end-2),2)+sum(InterferencePowerInBand_OFDM_15kHz(:,3:end-2),2)));
SIR_dB_Frequency_FBMC_15kHz(:,i_GuardBand)      = 10*log10(sum(SignalPowerInBand_FBMC_15kHz(:,5:end-4),2)./(sum(InterferencePower_120kHzTo15kHz_FBMC(:,5:end-4),2)+sum(InterferencePowerInBand_FBMC_15kHz(:,5:end-4),2)));
SIR_dB_Frequency_WOLA_15kHz(:,i_GuardBand)      = 10*log10(sum(SignalPowerInBand_WOLA_15kHz(:,3:end-2),2)./(sum(InterferencePower_120kHzTo15kHz_WOLA(:,3:end-2),2)+sum(InterferencePowerInBand_WOLA_15kHz(:,3:end-2),2)));
SIR_dB_Frequency_FOFDM_15kHz(:,i_GuardBand)     = 10*log10(sum(SignalPowerInBand_FOFDM_15kHz(:,3:end-2),2)./(sum(InterferencePower_120kHzTo15kHz_FOFDM(:,3:end-2),2)+sum(InterferencePowerInBand_FOFDM_15kHz(:,3:end-2),2)));
SIR_dB_Frequency_UFMC_15kHz(:,i_GuardBand)      = 10*log10(sum(SignalPowerInBand_UFMC_15kHz(:,3:end-2),2)./(sum(InterferencePower_120kHzTo15kHz_UFMC(:,3:end-2),2)+sum(InterferencePowerInBand_UFMC_15kHz(:,3:end-2),2)));

SIR_dB_Frequency_OFDM_120kHz(:,i_GuardBand)     = 10*log10(sum(SignalPowerInBand_OFDM_120kHz(:,3:end-2),2)./(sum(InterferencePower_15kHzTo120kHz_OFDM(:,3:end-2),2)+sum(InterferencePowerInBand_OFDM_120kHz(:,3:end-2),2)));
SIR_dB_Frequency_FBMC_120kHz(:,i_GuardBand)     = 10*log10(sum(SignalPowerInBand_FBMC_120kHz(:,5:end-4),2)./(sum(InterferencePower_15kHzTo120kHz_FBMC(:,5:end-4),2)+sum(InterferencePowerInBand_FBMC_120kHz(:,5:end-4),2)));
SIR_dB_Frequency_WOLA_120kHz(:,i_GuardBand)     = 10*log10(sum(SignalPowerInBand_WOLA_120kHz(:,3:end-2),2)./(sum(InterferencePower_15kHzTo120kHz_WOLA(:,3:end-2),2)+sum(InterferencePowerInBand_WOLA_120kHz(:,3:end-2),2)));
SIR_dB_Frequency_FOFDM_120kHz(:,i_GuardBand)    = 10*log10(sum(SignalPowerInBand_FOFDM_120kHz(:,3:end-2),2)./(sum(InterferencePower_15kHzTo120kHz_FOFDM(:,3:end-2),2)+sum(InterferencePowerInBand_FOFDM_120kHz(:,3:end-2),2)));
SIR_dB_Frequency_UFMC_120kHz(:,i_GuardBand)     = 10*log10(sum(SignalPowerInBand_UFMC_120kHz(:,3:end-2),2)./(sum(InterferencePower_15kHzTo120kHz_UFMC(:,3:end-2),2)+sum(InterferencePowerInBand_UFMC_120kHz(:,3:end-2),2)));

SIR_dB_Frequency_FBMC_15kHzFrom480kHz(:,i_GuardBand) = 10*log10(sum(SignalPowerInBand_FBMC_15kHz(:,5:end-4),2)./(sum(InterferencePower_480kHzTo15kHz_FBMC(:,5:end-4),2)+sum(InterferencePowerInBand_FBMC_15kHz(:,5:end-4),2)));
SIR_dB_Frequency_FBMC_480kHzFrom15kHz(:,i_GuardBand) = 10*log10(sum(SignalPowerInBand_FBMC_480kHz(:,5:end-4),2)./(sum(InterferencePower_15kHzTo480kHz_FBMC(:,5:end-4),2)+sum(InterferencePowerInBand_FBMC_480kHz(:,5:end-4),2)));


% Total Signal to Interference Ratio (time and frequency is averaged out and we remove the first 2 and last 2 symbols)
SIR_dB_Total_OFDM_15kHz(i_GuardBand)    = 10*log10(sum(sum(SignalPowerInBand_OFDM_15kHz(:,3:end-2),2))./(sum(sum(InterferencePower_120kHzTo15kHz_OFDM(:,3:end-2),2))+sum(sum(InterferencePowerInBand_OFDM_15kHz(:,3:end-2),2))));
SIR_dB_Total_FBMC_15kHz(i_GuardBand)    = 10*log10(sum(sum(SignalPowerInBand_FBMC_15kHz(:,5:end-4),2))./(sum(sum(InterferencePower_120kHzTo15kHz_FBMC(:,5:end-4),2))+sum(sum(InterferencePowerInBand_FBMC_15kHz(:,5:end-4),2))));
SIR_dB_Total_WOLA_15kHz(i_GuardBand)    = 10*log10(sum(sum(SignalPowerInBand_WOLA_15kHz(:,3:end-2),2))./(sum(sum(InterferencePower_120kHzTo15kHz_WOLA(:,3:end-2),2))+sum(sum(InterferencePowerInBand_WOLA_15kHz(:,3:end-2),2))));
SIR_dB_Total_FOFDM_15kHz(i_GuardBand)   = 10*log10(sum(sum(SignalPowerInBand_FOFDM_15kHz(:,3:end-2),2))./(sum(sum(InterferencePower_120kHzTo15kHz_FOFDM(:,3:end-2),2))+sum(sum(InterferencePowerInBand_FOFDM_15kHz(:,3:end-2),2))));
SIR_dB_Total_UFMC_15kHz(i_GuardBand)    = 10*log10(sum(sum(SignalPowerInBand_UFMC_15kHz(:,3:end-2),2))./(sum(sum(InterferencePower_120kHzTo15kHz_UFMC(:,3:end-2),2))+sum(sum(InterferencePowerInBand_UFMC_15kHz(:,3:end-2),2))));

SIR_dB_Total_OFDM_120kHz(i_GuardBand)   = 10*log10(sum(sum(SignalPowerInBand_OFDM_120kHz(:,3:end-2),2))./(sum(sum(InterferencePower_15kHzTo120kHz_OFDM(:,3:end-2),2))+sum(sum(InterferencePowerInBand_OFDM_120kHz(:,3:end-2),2))));
SIR_dB_Total_FBMC_120kHz(i_GuardBand)   = 10*log10(sum(sum(SignalPowerInBand_FBMC_120kHz(:,5:end-4),2))./(sum(sum(InterferencePower_15kHzTo120kHz_FBMC(:,5:end-4),2))+sum(sum(InterferencePowerInBand_FBMC_120kHz(:,5:end-4),2))));
SIR_dB_Total_WOLA_120kHz(i_GuardBand)   = 10*log10(sum(sum(SignalPowerInBand_WOLA_120kHz(:,3:end-2),2))./(sum(sum(InterferencePower_15kHzTo120kHz_WOLA(:,3:end-2),2))+sum(sum(InterferencePowerInBand_WOLA_120kHz(:,3:end-2),2))));
SIR_dB_Total_FOFDM_120kHz(i_GuardBand)  = 10*log10(sum(sum(SignalPowerInBand_FOFDM_120kHz(:,3:end-2),2))./(sum(sum(InterferencePower_15kHzTo120kHz_FOFDM(:,3:end-2),2))+sum(sum(InterferencePowerInBand_FOFDM_120kHz(:,3:end-2),2))));
SIR_dB_Total_UFMC_120kHz(i_GuardBand)   = 10*log10(sum(sum(SignalPowerInBand_UFMC_120kHz(:,3:end-2),2))./(sum(sum(InterferencePower_15kHzTo120kHz_UFMC(:,3:end-2),2))+sum(sum(InterferencePowerInBand_UFMC_120kHz(:,3:end-2),2))));

SIR_dB_Total_FBMC_15kHzFrom480kHz(i_GuardBand) = 10*log10(sum(sum(SignalPowerInBand_FBMC_15kHz(:,5:end-4),2))./(sum(sum(InterferencePower_480kHzTo15kHz_FBMC(:,5:end-4),2))+sum(sum(InterferencePowerInBand_FBMC_15kHz(:,5:end-4),2))));
SIR_dB_Total_FBMC_480kHzFrom15kHz(i_GuardBand) = 10*log10(sum(sum(SignalPowerInBand_FBMC_480kHz(:,5:end-4),2))./(sum(sum(InterferencePower_15kHzTo480kHz_FBMC(:,5:end-4),2))+sum(sum(InterferencePowerInBand_FBMC_480kHz(:,5:end-4),2))));

% Total Interference Both Users
SIR_dB_Total_OFDM(i_GuardBand) = 10*log10(...
    (sum(sum(SignalPowerInBand_OFDM_15kHz(:,3:end-2),2)) +  sum(sum(SignalPowerInBand_OFDM_120kHz(:,3:end-2),2))) ...
    ./( sum(sum(InterferencePower_120kHzTo15kHz_OFDM(:,3:end-2),2)) + sum(sum(InterferencePower_15kHzTo120kHz_OFDM(:,3:end-2),2)) ...
        +sum(sum(InterferencePowerInBand_OFDM_15kHz(:,3:end-2),2))  + sum(sum(InterferencePowerInBand_OFDM_120kHz(:,3:end-2),2)) ));
SIR_dB_Total_FBMC(i_GuardBand) = 10*log10(...
    (sum(sum(SignalPowerInBand_FBMC_15kHz(:,3:end-2),2)) +  sum(sum(SignalPowerInBand_FBMC_120kHz(:,3:end-2),2))) ...
    ./( sum(sum(InterferencePower_120kHzTo15kHz_FBMC(:,3:end-2),2)) + sum(sum(InterferencePower_15kHzTo120kHz_FBMC(:,3:end-2),2)) ...
        +sum(sum(InterferencePowerInBand_FBMC_15kHz(:,3:end-2),2))  + sum(sum(InterferencePowerInBand_FBMC_120kHz(:,3:end-2),2)) ));   
SIR_dB_Total_WOLA(i_GuardBand) = 10*log10(...
    (sum(sum(SignalPowerInBand_WOLA_15kHz(:,3:end-2),2)) +  sum(sum(SignalPowerInBand_WOLA_120kHz(:,3:end-2),2))) ...
    ./( sum(sum(InterferencePower_120kHzTo15kHz_WOLA(:,3:end-2),2)) + sum(sum(InterferencePower_15kHzTo120kHz_WOLA(:,3:end-2),2)) ...
        +sum(sum(InterferencePowerInBand_WOLA_15kHz(:,3:end-2),2))  + sum(sum(InterferencePowerInBand_WOLA_120kHz(:,3:end-2),2)) ));
SIR_dB_Total_FOFDM(i_GuardBand) = 10*log10(...
    (sum(sum(SignalPowerInBand_FOFDM_15kHz(:,3:end-2),2)) +  sum(sum(SignalPowerInBand_FOFDM_120kHz(:,3:end-2),2))) ...
    ./( sum(sum(InterferencePower_120kHzTo15kHz_FOFDM(:,3:end-2),2)) + sum(sum(InterferencePower_15kHzTo120kHz_FOFDM(:,3:end-2),2)) ...
        +sum(sum(InterferencePowerInBand_FOFDM_15kHz(:,3:end-2),2))  + sum(sum(InterferencePowerInBand_FOFDM_120kHz(:,3:end-2),2)) ));
SIR_dB_Total_UFMC(i_GuardBand) = 10*log10(...
    (sum(sum(SignalPowerInBand_UFMC_15kHz(:,3:end-2),2)) +  sum(sum(SignalPowerInBand_UFMC_120kHz(:,3:end-2),2))) ...
    ./( sum(sum(InterferencePower_120kHzTo15kHz_UFMC(:,3:end-2),2)) + sum(sum(InterferencePower_15kHzTo120kHz_UFMC(:,3:end-2),2)) ...
        +sum(sum(InterferencePowerInBand_UFMC_15kHz(:,3:end-2),2))  + sum(sum(InterferencePowerInBand_UFMC_120kHz(:,3:end-2),2)) ));
SIR_dB_Total_FBMC_480kHz(i_GuardBand) = 10*log10(...
    (sum(sum(SignalPowerInBand_FBMC_15kHz(:,3:end-2),2)) +  sum(sum(SignalPowerInBand_FBMC_480kHz(:,3:end-2),2))) ...
    ./( sum(sum(InterferencePower_480kHzTo15kHz_FBMC(:,3:end-2),2)) + sum(sum(InterferencePower_15kHzTo480kHz_FBMC(:,3:end-2),2)) ...
        +sum(sum(InterferencePowerInBand_FBMC_15kHz(:,3:end-2),2))  + sum(sum(InterferencePowerInBand_FBMC_480kHz(:,3:end-2),2)) ));   

 
%% The SIR according to the equation in the paper. Here we ignore the interband interference to keep the equation short
SIR_dB_PaperNotation_OFDM(i_GuardBand) =  10*log10((L_15kHz*(K_15kHz-4)+L_120kHz*(K_120kHz-4))/(norm(G_RX_OFDM_15kHz((2*L_15kHz+1):(end-2*L_15kHz),:)*G_TX_OFDM_120kHz_FrequencyShifted,'fro')^2+norm(G_RX_OFDM_120kHz_FrequencyShifted((2*L_120kHz+1):(end-2*L_120kHz),:)*G_TX_OFDM_15kHz,'fro')^2));
SIR_dB_PaperNotation_WOLA(i_GuardBand) =  10*log10((L_15kHz*(K_15kHz-4)+L_120kHz*(K_120kHz-4))/(norm(G_RX_WOLA_15kHz((2*L_15kHz+1):(end-2*L_15kHz),:)*G_TX_WOLA_120kHz_FrequencyShifted,'fro')^2+norm(G_RX_WOLA_120kHz_FrequencyShifted((2*L_120kHz+1):(end-2*L_120kHz),:)*G_TX_WOLA_15kHz,'fro')^2));
SIR_dB_PaperNotation_FOFDM(i_GuardBand) =  10*log10((L_15kHz*(K_15kHz-4)+L_120kHz*(K_120kHz-4))/(norm(G_RX_FOFDM_15kHz((2*L_15kHz+1):(end-2*L_15kHz),:)*G_TX_FOFDM_120kHz_FrequencyShifted,'fro')^2+norm(G_RX_FOFDM_120kHz_FrequencyShifted((2*L_120kHz+1):(end-2*L_120kHz),:)*G_TX_FOFDM_15kHz,'fro')^2));
SIR_dB_PaperNotation_UFMC(i_GuardBand) =  10*log10((L_15kHz*(K_15kHz-4)+L_120kHz*(K_120kHz-4))/(norm(G_RX_UFMC_15kHz((2*L_15kHz+1):(end-2*L_15kHz),:)*G_TX_UFMC_120kHz_FrequencyShifted,'fro')^2+norm(G_RX_UFMC_120kHz_FrequencyShifted((2*L_120kHz+1):(end-2*L_120kHz),:)*G_TX_UFMC_15kHz,'fro')^2));
SIR_dB_PaperNotation_FBMC(i_GuardBand) =  10*log10((L_15kHz*(2*K_15kHz-4)+L_120kHz*(K_120kHz*2-4))/(norm(real(G_RX_FBMC_15kHz((2*L_15kHz+1):(end-2*L_15kHz),:)*G_TX_FBMC_120kHz_FrequencyShifted),'fro')^2+norm(real(G_RX_FBMC_120kHz_FrequencyShifted((2*L_120kHz+1):(end-2*L_120kHz),:)*G_TX_FBMC_15kHz),'fro')^2));


    
% Calculate PSD
if M_NormalizedGuardBand(i_GuardBand)==NormalizedGuardBandForWhichPSDisPlotted
    f_OFDM  = (0:size(G_TX_OFDM_15kHz,1)-1)*SamplingRate/size(G_TX_OFDM_15kHz,1);
    f_FBMC  = (0:size(G_TX_FBMC_15kHz,1)-1)*SamplingRate/size(G_TX_FBMC_15kHz,1);
    f_WOLA  = (0:size(G_TX_WOLA_15kHz,1)-1)*SamplingRate/size(G_TX_WOLA_15kHz,1);
    f_FOFDM = (0:size(G_TX_FOFDM_15kHz,1)-1)*SamplingRate/size(G_TX_FOFDM_15kHz,1);
    f_UFMC  = (0:size(G_TX_UFMC_15kHz,1)-1)*SamplingRate/size(G_TX_UFMC_15kHz,1);
    
    PSD_OFDM_15kHz                      = sum(abs(fft(full(G_TX_OFDM_15kHz))).^2,2);
    PSD_OFDM_120kHz_FrequencyShifted    = sum(abs(fft(full(G_TX_OFDM_120kHz_FrequencyShifted))).^2,2);
    
    PSD_FBMC_15kHz                      = sum(abs(fft(full(G_TX_FBMC_15kHz))).^2,2);
    PSD_FBMC_120kHz_FrequencyShifted    = sum(abs(fft(full(G_TX_FBMC_120kHz_FrequencyShifted))).^2,2);    
    PSD_FBMC_480kHz_FrequencyShifted    = sum(abs(fft(full(G_TX_FBMC_480kHz_FrequencyShifted))).^2,2);    
   
    PSD_WOLA_15kHz                      = sum(abs(fft(full(G_TX_WOLA_15kHz))).^2,2);
    PSD_WOLA_120kHz_FrequencyShifted    = sum(abs(fft(full(G_TX_WOLA_120kHz_FrequencyShifted))).^2,2);      
    
    PSD_FOFDM_15kHz                     = sum(abs(fft(full(G_TX_FOFDM_15kHz))).^2,2);
    PSD_FOFDM_120kHz_FrequencyShifted   = sum(abs(fft(full(G_TX_FOFDM_120kHz_FrequencyShifted))).^2,2);      
        
    PSD_UFMC_15kHz                      = sum(abs(fft(full(G_TX_UFMC_15kHz))).^2,2);
    PSD_UFMC_120kHz_FrequencyShifted    = sum(abs(fft(full(G_TX_UFMC_120kHz_FrequencyShifted))).^2,2);        
end
    
disp([int2str(i_GuardBand/length(M_NormalizedGuardBand)*100) '%']);
end 


PSD_OFDM_15kHz_Normalized   = PSD_OFDM_15kHz/(sum(PSD_OFDM_15kHz)*SamplingRate/length(PSD_OFDM_15kHz));
PSD_OFDM_120kHz_Normalized  = PSD_OFDM_120kHz_FrequencyShifted/(sum(PSD_OFDM_120kHz_FrequencyShifted)*SamplingRate/length(PSD_OFDM_120kHz_FrequencyShifted));

PSD_FBMC_15kHz_Normalized   = PSD_FBMC_15kHz/(sum(PSD_FBMC_15kHz)*SamplingRate/length(PSD_FBMC_15kHz));
PSD_FBMC_120kHz_Normalized  = PSD_FBMC_120kHz_FrequencyShifted/(sum(PSD_FBMC_120kHz_FrequencyShifted)*SamplingRate/length(PSD_FBMC_120kHz_FrequencyShifted));
PSD_FBMC_480kHz_Normalized  = PSD_FBMC_480kHz_FrequencyShifted/(sum(PSD_FBMC_480kHz_FrequencyShifted)*SamplingRate/length(PSD_FBMC_480kHz_FrequencyShifted));

PSD_WOLA_15kHz_Normalized   = PSD_WOLA_15kHz/(sum(PSD_WOLA_15kHz)*SamplingRate/length(PSD_WOLA_15kHz));
PSD_WOLA_120kHz_Normalized  = PSD_WOLA_120kHz_FrequencyShifted/(sum(PSD_WOLA_120kHz_FrequencyShifted)*SamplingRate/length(PSD_WOLA_120kHz_FrequencyShifted));

PSD_FOFDM_15kHz_Normalized  = PSD_FOFDM_15kHz/(sum(PSD_FOFDM_15kHz)*SamplingRate/length(PSD_FOFDM_15kHz));
PSD_FOFDM_120kHz_Normalized = PSD_FOFDM_120kHz_FrequencyShifted/(sum(PSD_FOFDM_120kHz_FrequencyShifted)*SamplingRate/length(PSD_FOFDM_120kHz_FrequencyShifted));

PSD_UFMC_15kHz_Normalized   = PSD_UFMC_15kHz/(sum(PSD_UFMC_15kHz)*SamplingRate/length(PSD_UFMC_15kHz));
PSD_UFMC_120kHz_Normalized  = PSD_UFMC_120kHz_FrequencyShifted/(sum(PSD_UFMC_120kHz_FrequencyShifted)*SamplingRate/length(PSD_UFMC_120kHz_FrequencyShifted));


NormalizationFactor         = max(PSD_FBMC_15kHz_Normalized);
PSD_OFDM_15kHz_Normalized   = PSD_OFDM_15kHz_Normalized/NormalizationFactor;
PSD_OFDM_120kHz_Normalized  = PSD_OFDM_120kHz_Normalized/NormalizationFactor;

PSD_FBMC_15kHz_Normalized   = PSD_FBMC_15kHz_Normalized/NormalizationFactor;
PSD_FBMC_120kHz_Normalized  = PSD_FBMC_120kHz_Normalized/NormalizationFactor;
PSD_FBMC_480kHz_Normalized  = PSD_FBMC_480kHz_Normalized/NormalizationFactor;

PSD_WOLA_15kHz_Normalized   = PSD_WOLA_15kHz_Normalized/NormalizationFactor;
PSD_WOLA_120kHz_Normalized  = PSD_WOLA_120kHz_Normalized/NormalizationFactor;

PSD_FOFDM_15kHz_Normalized  = PSD_FOFDM_15kHz_Normalized/NormalizationFactor;
PSD_FOFDM_120kHz_Normalized = PSD_FOFDM_120kHz_Normalized/NormalizationFactor;

PSD_UFMC_15kHz_Normalized   = PSD_UFMC_15kHz_Normalized/NormalizationFactor;
PSD_UFMC_120kHz_Normalized  = PSD_UFMC_120kHz_Normalized/NormalizationFactor;


% Correct the negative frequencies in the PSD
N_OFDM  = length(f_OFDM);
N_FBMC  = length(f_FBMC);
N_WOLA  = length(f_WOLA);
N_FOFDM = length(f_FOFDM);
N_UFMC  = length(f_UFMC);

f_OFDM  = [f_OFDM(1:floor(N_OFDM/2)) -f_OFDM(ceil(N_OFDM/2)+1:-1:2)];
f_FBMC  = [f_FBMC(1:floor(N_FBMC/2)) -f_FBMC(ceil(N_FBMC/2)+1:-1:2)];
f_WOLA  = [f_WOLA(1:floor(N_WOLA/2)) -f_WOLA(ceil(N_WOLA/2)+1:-1:2)];
f_FOFDM = [f_FOFDM(1:floor(N_FOFDM/2)) -f_WOLA(ceil(N_FOFDM/2)+1:-1:2)];
f_UFMC  = [f_UFMC(1:floor(N_UFMC/2)) -f_UFMC(ceil(N_UFMC/2)+1:-1:2)];

f_OFDM  = circshift(f_OFDM,[0 ceil(N_OFDM/2)]);
f_FBMC  = circshift(f_FBMC,[0 ceil(N_FBMC/2)]);
f_WOLA  = circshift(f_WOLA,[0 ceil(N_WOLA/2)]);
f_FOFDM = circshift(f_FOFDM,[0 ceil(N_FOFDM/2)]);
f_UFMC  = circshift(f_UFMC,[0 ceil(N_UFMC/2)]);

PSD_OFDM_15kHz_Normalized   = circshift(PSD_OFDM_15kHz_Normalized,[ceil(N_OFDM/2) 0]);
PSD_OFDM_120kHz_Normalized  = circshift(PSD_OFDM_120kHz_Normalized,[ceil(N_OFDM/2) 0]);

PSD_FBMC_15kHz_Normalized   = circshift(PSD_FBMC_15kHz_Normalized,[ceil(N_FBMC/2) 0]);
PSD_FBMC_120kHz_Normalized  = circshift(PSD_FBMC_120kHz_Normalized,[ceil(N_FBMC/2) 0]);
PSD_FBMC_480kHz_Normalized  = circshift(PSD_FBMC_480kHz_Normalized,[ceil(N_FBMC/2) 0]);

PSD_WOLA_15kHz_Normalized   = circshift(PSD_WOLA_15kHz_Normalized,[ceil(N_WOLA/2) 0]);
PSD_WOLA_120kHz_Normalized  = circshift(PSD_WOLA_120kHz_Normalized,[ceil(N_WOLA/2) 0]);

PSD_FOFDM_15kHz_Normalized  = circshift(PSD_FOFDM_15kHz_Normalized,[ceil(N_FOFDM/2) 0]);
PSD_FOFDM_120kHz_Normalized = circshift(PSD_FOFDM_120kHz_Normalized,[ceil(N_FOFDM/2) 0]);

PSD_UFMC_15kHz_Normalized   = circshift(PSD_UFMC_15kHz_Normalized,[ceil(N_UFMC/2) 0]);
PSD_UFMC_120kHz_Normalized  = circshift(PSD_UFMC_120kHz_Normalized,[ceil(N_UFMC/2) 0]);



%% Plot Everthing
figure(1);
plot(M_NormalizedGuardBand,SIR_dB_Total_OFDM,'Color',ColorOFDM);
hold on;
plot(M_NormalizedGuardBand,SIR_dB_Total_FBMC,'Color',ColorFBMC);
plot(M_NormalizedGuardBand,SIR_dB_Total_FBMC_480kHz,':','Color',ColorFBMC);
plot(M_NormalizedGuardBand,SIR_dB_Total_WOLA,'Color',ColorWOLA);
plot(M_NormalizedGuardBand,SIR_dB_Total_FOFDM,'Color',ColorFOFDM);
plot(M_NormalizedGuardBand,SIR_dB_Total_UFMC,'Color',ColorUFMC);
ylim([0 60]);
xlabel('Normalized Guard Band, F_G/FL');
ylabel('Signal-to-Interference Ratio [dB]');

figure(11);
plot(M_NormalizedGuardBand,[SIR_dB_Total_OFDM.' SIR_dB_PaperNotation_OFDM'],'Color',ColorOFDM);
hold on;
plot(M_NormalizedGuardBand,[SIR_dB_Total_FBMC.' SIR_dB_PaperNotation_FBMC'],'Color',ColorFBMC);
plot(M_NormalizedGuardBand,SIR_dB_Total_FBMC_480kHz,':','Color',ColorFBMC);
plot(M_NormalizedGuardBand,[SIR_dB_Total_WOLA.' SIR_dB_PaperNotation_WOLA.'],'Color',ColorWOLA);
plot(M_NormalizedGuardBand,[SIR_dB_Total_FOFDM.' SIR_dB_PaperNotation_FOFDM.'],'Color',ColorFOFDM);
plot(M_NormalizedGuardBand,[SIR_dB_Total_UFMC.' SIR_dB_PaperNotation_UFMC.'],'Color',ColorUFMC);
ylim([0 60]);
xlabel('Normalized Guard Band, F_G/FL');
ylabel('Signal-to-Interference Ratio [dB]');
title('With and Without Inter-Band-Interference');


LF = L_15kHz*SubcarrierSpacing_15kHz;
figure(2);
plot(f_OFDM/LF,10*log10([PSD_OFDM_15kHz_Normalized PSD_OFDM_120kHz_Normalized]),'Color',ColorOFDM);hold on;
ylim([-60 2]);
xlim([-0.2 NormalizedGuardBandForWhichPSDisPlotted+2+0.2]);
plot([0 0], [-150,20],'color', [0.8 0.8 0.8]);
plot([1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[2 2], [-150,20],'color', [0.8 0.8 0.8]);
xlabel('Normalized Frequency f/FL');
ylabel('Power Spectral Density');
title('OFDM');

figure(3);
plot(f_FBMC/LF,10*log10([PSD_FBMC_15kHz_Normalized PSD_FBMC_120kHz_Normalized]),'Color',ColorFBMC);hold on;
ylim([-60 2]);
xlim([-0.2 NormalizedGuardBandForWhichPSDisPlotted+2+0.2]);
plot([0 0], [-150,20],'color', [0.8 0.8 0.8]);
plot([1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[2 2], [-150,20],'color', [0.8 0.8 0.8]);
xlabel('Normalized Frequency f/FL');
ylabel('Power Spectral Density');
title('FBMC 15kHz and 120kHz');

figure(4);
plot(f_FBMC/LF,10*log10([PSD_FBMC_15kHz_Normalized PSD_FBMC_480kHz_Normalized]),'Color',ColorFBMC);hold on;
ylim([-60 2]);
xlim([-0.2 NormalizedGuardBandForWhichPSDisPlotted+2+0.2]);
plot([0 0], [-150,20],'color', [0.8 0.8 0.8]);
plot([1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[2 2], [-150,20],'color', [0.8 0.8 0.8]);
xlabel('Normalized Frequency f/FL');
ylabel('Power Spectral Density');
title('FBMC 15kHz and 480kHz');

figure(5);
plot(f_WOLA/LF,10*log10([PSD_WOLA_15kHz_Normalized PSD_WOLA_120kHz_Normalized]),'Color',ColorWOLA);hold on;
ylim([-60 2]);
xlim([-0.2 NormalizedGuardBandForWhichPSDisPlotted+2+0.2]);
plot([0 0], [-150,20],'color', [0.8 0.8 0.8]);
plot([1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[2 2], [-150,20],'color', [0.8 0.8 0.8]);
xlabel('Normalized Frequency f/FL');
ylabel('Power Spectral Density');
title('WOLA');

figure(6);
plot(f_FOFDM/LF,10*log10([PSD_FOFDM_15kHz_Normalized PSD_FOFDM_120kHz_Normalized]),'Color',ColorFOFDM);hold on;
ylim([-60 2]);
xlim([-0.2 NormalizedGuardBandForWhichPSDisPlotted+2+0.2]);
plot([0 0], [-150,20],'color', [0.8 0.8 0.8]);
plot([1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[2 2], [-150,20],'color', [0.8 0.8 0.8]);
xlabel('Normalized Frequency f/FL');
ylabel('Power Spectral Density');
title('FOFDM');

figure(7);
plot(f_UFMC/LF,10*log10([PSD_UFMC_15kHz_Normalized PSD_UFMC_120kHz_Normalized]),'Color',ColorUFMC);hold on;
ylim([-60 2]);
xlim([-0.2 NormalizedGuardBandForWhichPSDisPlotted+2+0.2]);
plot([0 0], [-150,20],'color', [0.8 0.8 0.8]);
plot([1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[1 1], [-150,20],'color', [0.8 0.8 0.8]);
plot(NormalizedGuardBandForWhichPSDisPlotted+[2 2], [-150,20],'color', [0.8 0.8 0.8]);
xlabel('Normalized Frequency f/FL');
ylabel('Power Spectral Density');
title('UFMC');



