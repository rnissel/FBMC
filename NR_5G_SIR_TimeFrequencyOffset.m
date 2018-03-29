% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script calculates the Signal-to-Interference Ratio in case of a time
% and a frequency offset for different multicarrier modulation techniques:
%       1) OFDM with CP (worst spectral behaviour)
%       2) FBMC (best spectral properties, complex orthogonality is replaced by real orthogonality)
%       3) WOLA (windowed OFDM. The windowing is done at TX and RX)
%       4) FOFDM (filtered OFDM, sinc + Hann, filtering at TX and RX)
%       5) UFMC (filtered OFDM, subband-wise filtering, Dolph-Chebyshev window, cyclic prefix and zero padding are supported, filtering at TX and RX)
% Note that WOLA, FOFDM and UFMC have a high degree of freedom! Depending
% on the window/filter length and the CP/ZP length they perform differently.


clear; close all;

M_NormalizedTimeOffset              = 0:0.005:0.1;                      % Normalized time offset: t_off/T0=t_off*F
M_NormalizedFrequencyOffset         = 0:0.005:0.1;                      % Normalized frequency offset: f_off/F

ReceiveFilterForWOLA_UFMC_FOFDM     = true;                             % Set to false to investigate the effect if no receive filter is employed. Note that the TF-spacing could then be reducede a little bit (not done at the moment).
PrototypeFilterOQAM                 = 'PHYDYAS';                        % PHYDYAS, Hermite, RRC
OverlappingFactorFBMC               = 4;                                % Overlapping factor for FBMC. 
UFMC_ZeroPadding                    = true;                             % TRUE for Zero Padding (ZP) and FALSE for a conventional Cyclic Prefix (CP)

NrSubcarriers                       = 48;                               % Number of subcarrier                            
NrTimeSymbols                       = 8;                                % Number of complex valued symbols in time (FBMC: times two)
SubcarrierSpacing                   = 15e3;                             % Subcarrier spacing (Hz)

SamplingRate                        = SubcarrierSpacing*24*14;          % Samping rate (Samples/s). The factor 14 is helpful because then the CP fits perfectly (1/14=0.07 of the symbol duration)
dt                                  = 1/SamplingRate;

% Parameterset 1
% TF_FilteredAndWindowedOFDM = 1.09;                                      % Time frequency spacing for UFMC and FOFDM. Chosen so that the interference from filtering is small, that is, SIR>=65dB
% FilterLengthTXandRX_UFMC = 1/(14*SubcarrierSpacing);                    % 4.7탎, same as the CP in LTE. However, we need filtering at TX and RX => Total length of the filtering is two times larger
% FilterLengthTXandRX_FOFDM = 1.18*1/(14*SubcarrierSpacing);              % Increased by 1.2 because FOFDM is more robust than UFMC to ISI (due to the filter)

% Parameterset 2
% TF_FilteredAndWindowedOFDM = 1.28;                                    % Time frequency spacing for UFMC and FOFDM. Chosen so that the interference from filtering is small, that is, SIR>=65dB
% FilterLengthTXandRX_UFMC = 3/(14*SubcarrierSpacing);                  % 3 times as high as the CP in LTE to get a smaller OOB. 
% FilterLengthTXandRX_FOFDM = 1.33*3/(14*SubcarrierSpacing);            % Increased by 1.35 because FOFDM is more robust than UFMC to ISI (due to the filter)

% Parameterset 3 (perfect orthogonal)
% TF_FilteredAndWindowedOFDM = 1+2/(14);
% FilterLengthTXandRX_UFMC = 1/(14*SubcarrierSpacing);        
% FilterLengthTXandRX_FOFDM = 1/(14*SubcarrierSpacing);

% Parameterset 4 (LTE like, interference ~ 50)
TF_FilteredAndWindowedOFDM = 1+1/(14);
FilterLengthTXandRX_UFMC = 1/(14*SubcarrierSpacing);        
FilterLengthTXandRX_FOFDM = 1.1*1/(14*SubcarrierSpacing);

CP_LengthFilter = (TF_FilteredAndWindowedOFDM-1)*1/(SubcarrierSpacing);

ColorFBMC   = [0 0 1]*0.5;
ColorOFDM   = [1 0 0];
ColorWOLA   = [1 0 1];
ColorFOFDM  = [1 1 0]*0.7;
ColorUFMC   = [0 1 0]*0.5;


%% FBMC Object
FBMC = Modulation.FBMC(...
    NrSubcarriers,...                                                   % Number subcarriers
    NrTimeSymbols*2,...                                                 % Number FBMC symbols
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency first subcarrier (Hz)
    false,...                                                           % Transmit real valued signal
    [PrototypeFilterOQAM '-OQAM'],...                                   % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    OverlappingFactorFBMC, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                              % Initial phase shift
    true ...                                                            % Polyphase implementation
    );
FBMC_QAM = Modulation.FBMC(...
    NrSubcarriers,...                                                   % Number subcarriers
    NrTimeSymbols,...                                                   % Number FBMC symbols
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency first subcarrier (Hz)
    false,...                                                           % Transmit real valued signal
    'Hermite-QAM',...                                                   % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    OverlappingFactorFBMC/2, ...                                        % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                              % Initial phase shift
    true ...                                                            % Polyphase implementation
    );

%% Conventional OFDM Object
OFDM = Modulation.OFDM(...
    NrSubcarriers,...                                                   % Number subcarriers
    NrTimeSymbols,...                                                   % Number OFDM Symbols
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency first subcarrier (Hz)
    false,...                                                           % Transmit real valued signal
    1/(14*SubcarrierSpacing), ...                                       % Cyclic prefix length (s) 
    0 ...                                                               % Zero guard length (s)
    );

%% Filtered OFDM Object
FOFDM = Modulation.FOFDM(...
    NrSubcarriers,...                                                   % Number subcarriers
    NrTimeSymbols,...                                                   % Number FOFDM Symbols
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency first subcarrier (Hz)
    false,...                                                           % Transmit real valued signal
    0, ...                                                              % Cyclic prefix length (s) 
    0, ...                                                              % Zero guard length (s)
    FilterLengthTXandRX_FOFDM, ...                                      % Length of the transmit filter (s)
    ReceiveFilterForWOLA_UFMC_FOFDM*FilterLengthTXandRX_FOFDM, ...      % Length of the receive filter (s)
    CP_LengthFilter...                                                  % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
);

%% UFMC Oject
UFMC = Modulation.UFMC(...
    NrSubcarriers,...                                                   % Number subcarriers
    NrTimeSymbols,...                                                   % Number UFMC Symbols
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency first subcarrier (Hz)
    false,...                                                           % Transmit real valued signal
    0, ...                                                              % Cyclic prefix length (s) 
    0, ...                                                              % Zero guard length (s)
    FilterLengthTXandRX_UFMC, ...                                       % Length of the transmit filter (s)
    ReceiveFilterForWOLA_UFMC_FOFDM*FilterLengthTXandRX_UFMC, ...       % Length of the receive filter (s) 
    CP_LengthFilter, ...                                                % Length of the additional cyclic prefix (or zero guard symbol if ZP is used) in seconds (s). Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
    UFMC_ZeroPadding ...                                                % TRUE for Zero Padding (ZP) and FALSE for a conventional Cyclic Prefix (CP)
);

%% WOLA object
WOLA = Modulation.WOLA(...
    NrSubcarriers,...                                                   % Number subcarriers
    NrTimeSymbols,...                                                   % Number OFDM Symbols
    SubcarrierSpacing,...                                               % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency first subcarrier (Hz)
    false,...                                                           % Transmit real valued signal
    0, ...                                                              % Cyclic prefix length (s) 
    0, ...                                                              % Zero guard length (s)
    CP_LengthFilter/2, ...                                              % Length of the window overlapping (s) at the transmitter 
    ReceiveFilterForWOLA_UFMC_FOFDM*CP_LengthFilter/2 ...               % Length of the window overlapping (s) at the receiver
    );


%% Plot Information
fprintf('               |(complex)TF-Spacing| Bandwidth(LF)|   Time (KT)   \n');
fprintf('OFDM (with CP) |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', OFDM.PHY.TimeSpacing*OFDM.PHY.SubcarrierSpacing , OFDM.PHY.SubcarrierSpacing*OFDM.Nr.Subcarriers/1e6,OFDM.PHY.TimeSpacing*OFDM.Nr.MCSymbols/1e-6);
fprintf('WOLA           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', WOLA.PHY.TimeSpacing*WOLA.PHY.SubcarrierSpacing , WOLA.PHY.SubcarrierSpacing*WOLA.Nr.Subcarriers/1e6,WOLA.PHY.TimeSpacing*WOLA.Nr.MCSymbols/1e-6);
fprintf('FOFDM          |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', FOFDM.PHY.TimeSpacing*FOFDM.PHY.SubcarrierSpacing , FOFDM.PHY.SubcarrierSpacing*FOFDM.Nr.Subcarriers/1e6,FOFDM.PHY.TimeSpacing*FOFDM.Nr.MCSymbols/1e-6);
fprintf('UFMC           |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', UFMC.PHY.TimeSpacing*UFMC.PHY.SubcarrierSpacing , UFMC.PHY.SubcarrierSpacing*UFMC.Nr.Subcarriers/1e6,UFMC.PHY.TimeSpacing*UFMC.Nr.MCSymbols/1e-6);
fprintf('FBMC-OQAM      |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', FBMC.PHY.TimeSpacing*FBMC.PHY.SubcarrierSpacing*2 , FBMC.PHY.SubcarrierSpacing*FBMC.Nr.Subcarriers/1e6,FBMC.PHY.TimeSpacing*FBMC.Nr.MCSymbols/1e-6);
fprintf('FBMC-QAM       |%17.2f  |%8.2f MHz  | %8.2f 탎  | \n', FBMC_QAM.PHY.TimeSpacing*FBMC_QAM.PHY.SubcarrierSpacing , FBMC_QAM.PHY.SubcarrierSpacing*FBMC_QAM.Nr.Subcarriers/1e6,FBMC_QAM.PHY.TimeSpacing*FBMC_QAM.Nr.MCSymbols/1e-6);


%% Generate Transmit  Matrices
G_TX_OFDM       = sparse(OFDM.GetTXMatrix);
G_TX_FBMC       = sparse(FBMC.GetTXMatrix);
G_TX_WOLA       = sparse(WOLA.GetTXMatrix);
G_TX_FOFDM      = sparse(FOFDM.GetTXMatrix);
G_TX_UFMC       = sparse(UFMC.GetTXMatrix);
G_TX_FBMC_QAM   = sparse(FBMC_QAM.GetTXMatrix);

G_RX_OFDM       = sparse(OFDM.GetRXMatrix);
G_RX_FBMC       = sparse(FBMC.GetRXMatrix);
G_RX_WOLA       = sparse(WOLA.GetRXMatrix);
G_RX_FOFDM      = sparse(FOFDM.GetRXMatrix);
G_RX_UFMC       = sparse(UFMC.GetRXMatrix);
G_RX_FBMC_QAM   = sparse(FBMC_QAM.GetRXMatrix);


M_NormalizedTimeOffset_Samples      = round((M_NormalizedTimeOffset/SubcarrierSpacing)/dt);
M_NormalizedTimeOffset_Quantized    = (M_NormalizedTimeOffset_Samples*dt)*SubcarrierSpacing;
SIR_dB_TimeOffset_OFDM              = nan(length(M_NormalizedTimeOffset_Samples),1);
SIR_dB_TimeOffset_FBMC              = nan(length(M_NormalizedTimeOffset_Samples),1);
SIR_dB_TimeOffset_WOLA              = nan(length(M_NormalizedTimeOffset_Samples),1);
SIR_dB_TimeOffset_FOFDM             = nan(length(M_NormalizedTimeOffset_Samples),1);
SIR_dB_TimeOffset_UFMC              = nan(length(M_NormalizedTimeOffset_Samples),1);
SIR_dB_TimeOffset_FBMC_QAM          = nan(length(M_NormalizedTimeOffset_Samples),1);
disp('Time offset...');
for i_timeOffset = 1:length(M_NormalizedTimeOffset_Samples)
    timeOffset  = M_NormalizedTimeOffset_Samples(i_timeOffset);
    D_OFDM      = circshift(G_RX_OFDM     ,[0,-timeOffset]) * G_TX_OFDM;
    D_WOLA      = circshift(G_RX_WOLA     ,[0,-timeOffset]) * G_TX_WOLA;
    D_FOFDM     = circshift(G_RX_FOFDM    ,[0,-timeOffset]) * G_TX_FOFDM;
    D_UFMC      = circshift(G_RX_UFMC     ,[0,-timeOffset]) * G_TX_UFMC;
    D_FBMC_QAM  = circshift(G_RX_FBMC_QAM ,[0,-timeOffset]) * G_TX_FBMC_QAM;
    
    D_FBMC      = circshift(G_RX_FBMC     ,[0,-timeOffset]) * G_TX_FBMC;
    % Compensate the phase in FBMC!
    D_FBMC      = D_FBMC.*repmat(1./diag(D_FBMC).*abs(diag(D_FBMC)),[1,size(D_FBMC,1)]);
    
    SignalPowerInBand_OFDM              = full(reshape(abs(diag(D_OFDM)).^2,OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols));
    InterferencePowerInBand_OFDM        = reshape(sum(abs(D_OFDM).^2,2),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols) - SignalPowerInBand_OFDM;

    SignalPowerInBand_FBMC              = full(reshape(abs(diag(real(D_FBMC))).^2,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols));
    InterferencePowerInBand_FBMC        = reshape(sum(abs(real(D_FBMC)).^2,2),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols) - SignalPowerInBand_FBMC;

    SignalPowerInBand_WOLA              = full(reshape(abs(diag(D_WOLA)).^2,WOLA.Nr.Subcarriers,WOLA.Nr.MCSymbols));
    InterferencePowerInBand_WOLA        = reshape(sum(abs(D_WOLA).^2,2),WOLA.Nr.Subcarriers,WOLA.Nr.MCSymbols) - SignalPowerInBand_WOLA;

    SignalPowerInBand_FOFDM             = full(reshape(abs(diag(D_FOFDM)).^2,FOFDM.Nr.Subcarriers,FOFDM.Nr.MCSymbols));
    InterferencePowerInBand_FOFDM       = reshape(sum(abs(D_FOFDM).^2,2),FOFDM.Nr.Subcarriers,FOFDM.Nr.MCSymbols) - SignalPowerInBand_FOFDM;

    SignalPowerInBand_UFMC              = full(reshape(abs(diag(D_UFMC)).^2,UFMC.Nr.Subcarriers,UFMC.Nr.MCSymbols));
    InterferencePowerInBand_UFMC        = reshape(sum(abs(D_UFMC).^2,2),UFMC.Nr.Subcarriers,UFMC.Nr.MCSymbols) - SignalPowerInBand_UFMC;

    SignalPowerInBand_FBMC_QAM          = full(reshape(abs(diag(D_FBMC_QAM)).^2,FBMC_QAM.Nr.Subcarriers,FBMC_QAM.Nr.MCSymbols));
    InterferencePowerInBand_FBMC_QAM    = reshape(sum(abs(D_FBMC_QAM).^2,2),FBMC_QAM.Nr.Subcarriers,FBMC_QAM.Nr.MCSymbols) - SignalPowerInBand_FBMC_QAM;

    
    SIR_dB_TimeOffset_OFDM(i_timeOffset)        = 10*log10(sum(sum(SignalPowerInBand_OFDM(:,:)))./sum(sum(InterferencePowerInBand_OFDM(:,:))));
    SIR_dB_TimeOffset_FBMC(i_timeOffset)        = 10*log10(sum(sum(SignalPowerInBand_FBMC(:,:)))./sum(sum(InterferencePowerInBand_FBMC(:,:))));
    SIR_dB_TimeOffset_WOLA(i_timeOffset)        = 10*log10(sum(sum(SignalPowerInBand_WOLA(:,:)))./sum(sum(InterferencePowerInBand_WOLA(:,:))));
    SIR_dB_TimeOffset_FOFDM(i_timeOffset)       = 10*log10(sum(sum(SignalPowerInBand_FOFDM(:,:)))./sum(sum(InterferencePowerInBand_FOFDM(:,:))));
    SIR_dB_TimeOffset_UFMC(i_timeOffset)        = 10*log10(sum(sum(SignalPowerInBand_UFMC(:,:)))./sum(sum(InterferencePowerInBand_UFMC(:,:))));
    SIR_dB_TimeOffset_FBMC_QAM(i_timeOffset)    = 10*log10(sum(sum(SignalPowerInBand_FBMC_QAM(:,:)))./sum(sum(InterferencePowerInBand_FBMC_QAM(:,:))));
    
    disp([int2str(i_timeOffset/length(M_NormalizedTimeOffset_Samples)*100) '%']);
end
% Plot Stuff
figure();
plot(M_NormalizedTimeOffset_Quantized , SIR_dB_TimeOffset_OFDM     , 'Color' , ColorOFDM);hold on;
plot(M_NormalizedTimeOffset_Quantized , SIR_dB_TimeOffset_FBMC     , 'Color' , ColorFBMC);hold on;
plot(M_NormalizedTimeOffset_Quantized , SIR_dB_TimeOffset_WOLA     , 'Color' , ColorWOLA);hold on;
plot(M_NormalizedTimeOffset_Quantized , SIR_dB_TimeOffset_FOFDM    , 'Color' , ColorFOFDM);hold on;
plot(M_NormalizedTimeOffset_Quantized , SIR_dB_TimeOffset_UFMC     , 'Color' , ColorUFMC);hold on;
plot(M_NormalizedTimeOffset_Quantized , SIR_dB_TimeOffset_FBMC_QAM , 'blue') ; hold on;
ylim([0 50]);
xlim([M_NormalizedTimeOffset_Quantized(1) M_NormalizedTimeOffset_Quantized(end)]);
legend({'CP-OFDM','FBMC-OQAM','WOLA','FOFDM','UFMC','FBMC-QAM'});
xlabel('Normalized Time Offset, tF');
ylabel('Signal-to-Interference Ratio [dB]');


t_OFDM      =  (0:size(G_TX_OFDM,1)-1)*dt;
t_FBMC      =  (0:size(G_TX_FBMC,1)-1)*dt;
t_WOLA      =  (0:size(G_TX_WOLA,1)-1)*dt;
t_FOFDM     =  (0:size(G_TX_FOFDM,1)-1)*dt;
t_UFMC      =  (0:size(G_TX_UFMC,1)-1)*dt;
t_FBMC_QAM  =  (0:size(G_TX_FBMC_QAM,1)-1)*dt;


SIR_dB_FrequencyOffset_OFDM     = nan(length(M_NormalizedFrequencyOffset),1);
SIR_dB_FrequencyOffset_FBMC     = nan(length(M_NormalizedFrequencyOffset),1);
SIR_dB_FrequencyOffset_WOLA     = nan(length(M_NormalizedFrequencyOffset),1);
SIR_dB_FrequencyOffset_FOFDM    = nan(length(M_NormalizedFrequencyOffset),1);
SIR_dB_FrequencyOffset_UFMC     = nan(length(M_NormalizedFrequencyOffset),1);
SIR_dB_FrequencyOffset_FBMC_QAM = nan(length(M_NormalizedFrequencyOffset),1);
disp('Frequency offset...');
for i_frequencyOffset = 1:length(M_NormalizedFrequencyOffset)
    NormalizedFrequencyOffset = M_NormalizedFrequencyOffset(i_frequencyOffset);

    G_RX_OFDM_FrequencyOffset       = bsxfun(@times,G_RX_OFDM,exp(-1j*2*pi*(NormalizedFrequencyOffset*SubcarrierSpacing)*t_OFDM));
    G_RX_FBMC_FrequencyOffset       = bsxfun(@times,G_RX_FBMC,exp(-1j*2*pi*(NormalizedFrequencyOffset*SubcarrierSpacing)*t_FBMC));
    G_RX_WOLA_FrequencyOffset       = bsxfun(@times,G_RX_WOLA,exp(-1j*2*pi*(NormalizedFrequencyOffset*SubcarrierSpacing)*t_WOLA));
    G_RX_FOFDM_FrequencyOffset      = bsxfun(@times,G_RX_FOFDM,exp(-1j*2*pi*(NormalizedFrequencyOffset*SubcarrierSpacing)*t_FOFDM));
    G_RX_UFMC_FrequencyOffset       = bsxfun(@times,G_RX_UFMC,exp(-1j*2*pi*(NormalizedFrequencyOffset*SubcarrierSpacing)*t_UFMC));
    G_RX_FBMC_QAM_FrequencyOffset   = bsxfun(@times,G_RX_FBMC_QAM,exp(-1j*2*pi*(NormalizedFrequencyOffset*SubcarrierSpacing)*t_FBMC_QAM));


    D_OFDM          = (G_RX_OFDM_FrequencyOffset)*G_TX_OFDM;
    D_WOLA          = (G_RX_WOLA_FrequencyOffset)*G_TX_WOLA;
    D_FOFDM         = (G_RX_FOFDM_FrequencyOffset)*G_TX_FOFDM;
    D_UFMC          = (G_RX_UFMC_FrequencyOffset)*G_TX_UFMC;
    D_FBMC_QAM      = (G_RX_FBMC_QAM_FrequencyOffset)*G_TX_FBMC_QAM;
    
    D_FBMC          = (G_RX_FBMC_FrequencyOffset)*G_TX_FBMC;
    % Compensate the phase in FBMC!
    D_FBMC          = D_FBMC.*repmat(1./diag(D_FBMC).*abs(diag(D_FBMC)),[1,size(D_FBMC,1)]);

    SignalPowerInBand_OFDM              = full(reshape(abs(diag(D_OFDM)).^2,OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols));
    InterferencePowerInBand_OFDM        = reshape(sum(abs(D_OFDM).^2,2),OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols) - SignalPowerInBand_OFDM;

    SignalPowerInBand_FBMC              = full(reshape(abs(diag(real(D_FBMC))).^2,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols));
    InterferencePowerInBand_FBMC        = reshape(sum(abs(real(D_FBMC)).^2,2),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols) - SignalPowerInBand_FBMC;

    SignalPowerInBand_WOLA              = full(reshape(abs(diag(D_WOLA)).^2,WOLA.Nr.Subcarriers,WOLA.Nr.MCSymbols));
    InterferencePowerInBand_WOLA        = reshape(sum(abs(D_WOLA).^2,2),WOLA.Nr.Subcarriers,WOLA.Nr.MCSymbols) - SignalPowerInBand_WOLA;

    SignalPowerInBand_FOFDM             = full(reshape(abs(diag(D_FOFDM)).^2,FOFDM.Nr.Subcarriers,FOFDM.Nr.MCSymbols));
    InterferencePowerInBand_FOFDM       = reshape(sum(abs(D_FOFDM).^2,2),FOFDM.Nr.Subcarriers,FOFDM.Nr.MCSymbols) - SignalPowerInBand_FOFDM;

    SignalPowerInBand_UFMC              = full(reshape(abs(diag(D_UFMC)).^2,UFMC.Nr.Subcarriers,UFMC.Nr.MCSymbols));
    InterferencePowerInBand_UFMC        = reshape(sum(abs(D_UFMC).^2,2),UFMC.Nr.Subcarriers,UFMC.Nr.MCSymbols) - SignalPowerInBand_UFMC;

    SignalPowerInBand_FBMC_QAM          = full(reshape(abs(diag(D_FBMC_QAM)).^2,FBMC_QAM.Nr.Subcarriers,FBMC_QAM.Nr.MCSymbols));
    InterferencePowerInBand_FBMC_QAM    = reshape(sum(abs(D_FBMC_QAM).^2,2),FBMC_QAM.Nr.Subcarriers,FBMC_QAM.Nr.MCSymbols) - SignalPowerInBand_FBMC_QAM;
    
    SIR_dB_FrequencyOffset_OFDM(i_frequencyOffset)      = 10*log10(sum(sum(SignalPowerInBand_OFDM(:,:)))./sum(sum(InterferencePowerInBand_OFDM(:,:))));
    SIR_dB_FrequencyOffset_FBMC(i_frequencyOffset)      = 10*log10(sum(sum(SignalPowerInBand_FBMC(:,:)))./sum(sum(InterferencePowerInBand_FBMC(:,:))));
    SIR_dB_FrequencyOffset_WOLA(i_frequencyOffset)      = 10*log10(sum(sum(SignalPowerInBand_WOLA(:,:)))./sum(sum(InterferencePowerInBand_WOLA(:,:))));
    SIR_dB_FrequencyOffset_FOFDM(i_frequencyOffset)     = 10*log10(sum(sum(SignalPowerInBand_FOFDM(:,:)))./sum(sum(InterferencePowerInBand_FOFDM(:,:))));
    SIR_dB_FrequencyOffset_UFMC(i_frequencyOffset)      = 10*log10(sum(sum(SignalPowerInBand_UFMC(:,:)))./sum(sum(InterferencePowerInBand_UFMC(:,:))));
    SIR_dB_FrequencyOffset_FBMC_QAM(i_frequencyOffset)  = 10*log10(sum(sum(SignalPowerInBand_FBMC_QAM(:,:)))./sum(sum(InterferencePowerInBand_FBMC_QAM(:,:))));
    
    disp([int2str(i_frequencyOffset/length(M_NormalizedFrequencyOffset)*100) '%']);
end
% Plot Stuff
figure();
plot(M_NormalizedFrequencyOffset , SIR_dB_FrequencyOffset_OFDM     , 'Color' , ColorOFDM);hold on;
plot(M_NormalizedFrequencyOffset , SIR_dB_FrequencyOffset_FBMC     , 'Color' , ColorFBMC);hold on;
plot(M_NormalizedFrequencyOffset , SIR_dB_FrequencyOffset_WOLA     , 'Color' , ColorWOLA);hold on;
plot(M_NormalizedFrequencyOffset , SIR_dB_FrequencyOffset_FOFDM    , 'Color' , ColorFOFDM);hold on;
plot(M_NormalizedFrequencyOffset , SIR_dB_FrequencyOffset_UFMC     , 'Color' , ColorUFMC);hold on;
plot(M_NormalizedFrequencyOffset , SIR_dB_FrequencyOffset_FBMC_QAM , 'blue');hold on;
ylim([0 50]);
xlim([M_NormalizedFrequencyOffset(1) M_NormalizedFrequencyOffset(end)]);
legend({'CP-OFDM','FBMC-OQAM','WOLA','FOFDM','UFMC','FBMC-QAM'});
xlabel('Normalized Frequency Offset, f/F');
ylabel('Signal-to-Interference Ratio [dB]');


