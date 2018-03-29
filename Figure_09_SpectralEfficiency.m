% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script plots the time-frequency efficiency of FBMC and f-OFDM.
% See Figure 9 in the paper (slightly changed parameters to improve the
% simulation time)

clear; close all;

%% Parameters
M_L = [1:1:10 12:4:20 20:8:128];   % [1:1:125]                          % For loop: number of subcarriers

EnergyThreshold         = 0.9999;                                       % 0.9999/(1-0.9999)=40dB

K                       = 1;                                            % Number of complex symbols in time
SubcarrierSpacing       = 15e3;                                         % Subcarrier spacing (15kHz, same as LTE)
N_FFT                   = 1*1024;  % 10*1024                            % FFT size
OverlappingFactor       = 4;                                            % Overlapping factor in FBMC
AddedZeroFactorForPSD   = 2;       %10                                  % Increase the number of zeros in time to improve the resolution in frequency => smoother curves. The factor is a multiple T0=1/F.


FOFDM_TF1 = 1+1/14;                                                     % Time-Frequency spacing 1 for FOFDM
FOFDM_TF2 = 1+2/14;                                                     % Time-Frequency spacing 2 for FOFDM
FOFDM_TF3 = 1+3/14;                                                     % Time-Frequency spacing 3 for FOFDM


%% Start Calculation 
SamplingRate            = SubcarrierSpacing*N_FFT;                      % Sampling rate
dt                      = 1/SamplingRate;
for i_L = 1:length(M_L)

    L = M_L(i_L);                                                       % Number of subcarriers
    IntermediateSubcarrier = round(N_FFT/2-L/2);

    %% FBMC Object
    FBMC_Hermite = Modulation.FBMC(...
        L,...                                                           % Number subcarriers
        2*K,...                                                         % Number FBMC symbols
        SubcarrierSpacing,...                                           % Subcarrier spacing (Hz)
        SamplingRate,...                                                % Sampling rate (Samples/s)
        SubcarrierSpacing*IntermediateSubcarrier,...                    % Intermediate frequency first subcarrier (Hz)
        false,...                                                       % Transmit real valued signal
        'Hermite-OQAM',...                                              % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
        OverlappingFactor, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
        0, ...                                                          % Initial phase shift
        true ...                                                        % Polyphase implementation
        );
    FBMC_PHYDYAS = Modulation.FBMC(...
        L,...                                                           % Number subcarriers
        2*K,...                                                         % Number FBMC symbols
        SubcarrierSpacing,...                                           % Subcarrier spacing (Hz)
        SamplingRate,...                                                % Sampling rate (Samples/s)
        SubcarrierSpacing*IntermediateSubcarrier,...                    % Intermediate frequency first subcarrier (Hz)
        false,...                                                       % Transmit real valued signal
        'PHYDYAS-OQAM',...                                              % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
        OverlappingFactor, ...                                          % Overlapping factor (also determines oversampling in the frequency domain)
        0, ...                                                          % Initial phase shift
        true ...                                                        % Polyphase implementation
        );

    %% f-OFDM Object
    FOFDM1 = Modulation.FOFDM(...
        L,...                                                           % Number subcarriers
        K,...                                                           % Number FOFDM Symbols
        SubcarrierSpacing,...                                           % Subcarrier spacing (Hz)
        SamplingRate,...                                                % Sampling rate (Samples/s)
        SubcarrierSpacing*IntermediateSubcarrier,...                    % Intermediate frequency first subcarrier (Hz)
        false,...                                                       % Transmit real valued signal
        0, ...                                                          % Cyclic prefix length (s) 
        AddedZeroFactorForPSD*1/SubcarrierSpacing, ...                  % Zero guard length (s); Large zero guard to increase the resolution of the power spectral density
        (FOFDM_TF1-1)*1/SubcarrierSpacing, ...                          % Length of the transmit filter (s)
        (FOFDM_TF1-1)*1/SubcarrierSpacing, ...                          % Length of the receive filter (s) 
        (FOFDM_TF1-1)*1/SubcarrierSpacing...                            % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
    );
    FOFDM2 = Modulation.FOFDM(...
        L,...                                                           % Number subcarriers
        K,...                                                           % Number FOFDM Symbols
        SubcarrierSpacing,...                                           % Subcarrier spacing (Hz)
        SamplingRate,...                                                % Sampling rate (Samples/s)
        SubcarrierSpacing*IntermediateSubcarrier,...                    % Intermediate frequency first subcarrier (Hz)
        false,...                                                       % Transmit real valued signal
        0, ...                                                          % Cyclic prefix length (s) 
        AddedZeroFactorForPSD*1/SubcarrierSpacing, ...                  % Zero guard length (s); Large zero guard to increase the resolution of the power spectral density
        (FOFDM_TF2-1)*1/SubcarrierSpacing, ...                          % Length of the transmit filter (s)
        (FOFDM_TF2-1)*1/SubcarrierSpacing, ...                          % Length of the receive filter (s) 
        (FOFDM_TF2-1)/SubcarrierSpacing...                              % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
    );
    FOFDM3 = Modulation.FOFDM(...
        L,...                                                           % Number subcarriers
        K,...                                                           % Number FOFDM Symbols
        SubcarrierSpacing,...                                           % Subcarrier spacing (Hz)
        SamplingRate,...                                                % Sampling rate (Samples/s)
        SubcarrierSpacing*IntermediateSubcarrier,...                    % Intermediate frequency first subcarrier (Hz)
        false,...                                                       % Transmit real valued signal
        0, ...                                                          % Cyclic prefix length (s) 
        AddedZeroFactorForPSD*1/SubcarrierSpacing, ...                  % Zero guard length (s); Large zero guard to increase the resolution of the power spectral density
        (FOFDM_TF3-1)*1/SubcarrierSpacing, ...                          % Length of the transmit filter (s)
        (FOFDM_TF3-1)*1/SubcarrierSpacing, ...                          % Length of the receive filter (s) 
        (FOFDM_TF3-1)*1/SubcarrierSpacing...                            % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
    );

    % Calculate the expected power over time
    [PS_FBMC_Hermite,~] = FBMC_Hermite.PlotTransmitPower;
    [PS_FBMC_PHYDYAS,~] = FBMC_PHYDYAS.PlotTransmitPower;
    [PS_FOFDM1,~]       = FOFDM1.PlotTransmitPower;
    [PS_FOFDM2,~]       = FOFDM2.PlotTransmitPower;
    [PS_FOFDM3,~]       = FOFDM3.PlotTransmitPower;

    % Calculate the required time according to the threshold
    convPS_FBMC_Hermite                     = conv(PS_FBMC_Hermite,ones(size(PS_FBMC_Hermite)));
    convPS_FBMC_Hermite(ceil(end/2)+1:end)  = [];
    Threshold_FBMC_Hermite_Upper            = (convPS_FBMC_Hermite(end)-convPS_FBMC_Hermite(end)*(1-EnergyThreshold)/2);
    Threshold_FBMC_Hermite_Lower            = (convPS_FBMC_Hermite(end)*(1-EnergyThreshold)/2);
    RequiredTime_FBMC_Hermite               = sum((convPS_FBMC_Hermite<Threshold_FBMC_Hermite_Upper) & (convPS_FBMC_Hermite>Threshold_FBMC_Hermite_Lower))*(dt);

    convPS_FBMC_PHYDYAS                     = conv(PS_FBMC_PHYDYAS,ones(size(PS_FBMC_PHYDYAS)));
    convPS_FBMC_PHYDYAS(ceil(end/2)+1:end)  = [];
    Threshold_FBMC_PHYDYAS_Upper            = (convPS_FBMC_PHYDYAS(end)-convPS_FBMC_PHYDYAS(end)*(1-EnergyThreshold)/2);
    Threshold_FBMC_PHYDYAS_Lower            = (convPS_FBMC_PHYDYAS(end)*(1-EnergyThreshold)/2);
    RequiredTime_FBMC_PHYDYAS               = sum((convPS_FBMC_PHYDYAS<Threshold_FBMC_PHYDYAS_Upper) & (convPS_FBMC_PHYDYAS>Threshold_FBMC_PHYDYAS_Lower))*(dt);

    convPS_FOFDM1                     = conv(PS_FOFDM1,ones(size(PS_FOFDM1)));
    convPS_FOFDM1(ceil(end/2)+1:end)  = [];
    Threshold_FOFDM1_Upper            = (convPS_FOFDM1(end)-convPS_FOFDM1(end)*(1-EnergyThreshold)/2);
    Threshold_FOFDM1_Lower            = (convPS_FOFDM1(end)*(1-EnergyThreshold)/2);
    RequiredTime_FOFDM1               = sum((convPS_FOFDM1<Threshold_FOFDM1_Upper) & (convPS_FOFDM1>Threshold_FOFDM1_Lower))*(dt);

    convPS_FOFDM2                     = conv(PS_FOFDM2,ones(size(PS_FOFDM2)));
    convPS_FOFDM2(ceil(end/2)+1:end)  = [];
    Threshold_FOFDM2_Upper            = (convPS_FOFDM2(end)-convPS_FOFDM2(end)*(1-EnergyThreshold)/2);
    Threshold_FOFDM2_Lower            = (convPS_FOFDM2(end)*(1-EnergyThreshold)/2);
    RequiredTime_FOFDM2               = sum((convPS_FOFDM2<Threshold_FOFDM2_Upper) & (convPS_FOFDM2>Threshold_FOFDM2_Lower))*(dt);

    convPS_FOFDM3                     = conv(PS_FOFDM3,ones(size(PS_FOFDM3)));
    convPS_FOFDM3(ceil(end/2)+1:end)  = [];
    Threshold_FOFDM3_Upper            = (convPS_FOFDM3(end)-convPS_FOFDM3(end)*(1-EnergyThreshold)/2);
    Threshold_FOFDM3_Lower            = (convPS_FOFDM3(end)*(1-EnergyThreshold)/2);
    RequiredTime_FOFDM3               = sum((convPS_FOFDM3<Threshold_FOFDM3_Upper) & (convPS_FOFDM3>Threshold_FOFDM3_Lower))*(dt);


    % Calculate the guard time, T_G
    GuardTime_FBMC_Hermite(i_L) = (RequiredTime_FBMC_Hermite-FBMC_Hermite.PHY.TimeSpacing*FBMC_Hermite.Nr.MCSymbols);
    GuardTime_FBMC_PHYDYAS(i_L) = (RequiredTime_FBMC_PHYDYAS-FBMC_PHYDYAS.PHY.TimeSpacing*FBMC_PHYDYAS.Nr.MCSymbols);
    GuardTime_FOFDM1(i_L)       = (RequiredTime_FOFDM1-FOFDM1.PHY.TimeSpacing*FOFDM1.Nr.MCSymbols);
    GuardTime_FOFDM2(i_L)       = (RequiredTime_FOFDM2-FOFDM2.PHY.TimeSpacing*FOFDM2.Nr.MCSymbols);
    GuardTime_FOFDM3(i_L)       = (RequiredTime_FOFDM3-FOFDM3.PHY.TimeSpacing*FOFDM3.Nr.MCSymbols);


    % Calculate the power spectral density
    [PSD_FOFDM1,f_FOFDM1] = FOFDM1.PlotPowerSpectralDensityUncorrelatedData;
    [PSD_FOFDM2,f_FOFDM2] = FOFDM2.PlotPowerSpectralDensityUncorrelatedData;
    [PSD_FOFDM3,f_FOFDM3] = FOFDM3.PlotPowerSpectralDensityUncorrelatedData;

    % Calulate the PSD for FBMC where we add zeros to increase the frequency resolution
    ZerosToAdd = zeros(AddedZeroFactorForPSD*FBMC_Hermite.Implementation.FFTSize,1);
    PSD_FBMC_Hermite = zeros(FBMC_Hermite.Nr.SamplesTotal+length(ZerosToAdd),1);
    PSD_FBMC_PHYDYAS = zeros(FBMC_Hermite.Nr.SamplesTotal+length(ZerosToAdd),1);
    for i_l = 1:L
        V = zeros(L,2*K);
        V(i_l,K)=1;
        PSD_FBMC_Hermite = PSD_FBMC_Hermite+abs(fft([FBMC_Hermite.Modulation(V);ZerosToAdd])).^2;
        PSD_FBMC_PHYDYAS = PSD_FBMC_PHYDYAS+abs(fft([FBMC_PHYDYAS.Modulation(V);ZerosToAdd])).^2;
    end
    f_FBMC = (0:length(PSD_FBMC_Hermite)-1)*1/(length(PSD_FBMC_Hermite)*dt);

    % Calculate the required frequency resources according to the threshold
    convPSD_FBMC_Hermite=conv(PSD_FBMC_Hermite,ones(size(PSD_FBMC_Hermite)));
    convPSD_FBMC_Hermite(ceil(end/2)+1:end) = [];
    Threshold_FBMC_Hermite_Upper            = (convPSD_FBMC_Hermite(end)-convPSD_FBMC_Hermite(end)*(1-EnergyThreshold)/2);
    Threshold_FBMC_Hermite_Lower            = (convPSD_FBMC_Hermite(end)*(1-EnergyThreshold)/2);
    RequiredBandwidth_FBMC_Hermite = sum((convPSD_FBMC_Hermite<Threshold_FBMC_Hermite_Upper) & (convPSD_FBMC_Hermite>Threshold_FBMC_Hermite_Lower))*(f_FBMC(2)-f_FBMC(1));

    convPSD_FBMC_PHYDYAS=conv(PSD_FBMC_PHYDYAS,ones(size(PSD_FBMC_PHYDYAS)));
    convPSD_FBMC_PHYDYAS(ceil(end/2)+1:end) = [];
    Threshold_FBMC_PHYDYAS_Upper            = (convPSD_FBMC_PHYDYAS(end)-convPSD_FBMC_PHYDYAS(end)*(1-EnergyThreshold)/2);
    Threshold_FBMC_PHYDYAS_Lower            = (convPSD_FBMC_PHYDYAS(end)*(1-EnergyThreshold)/2);
    RequiredBandwidth_FBMC_PHYDYAS          = sum((convPSD_FBMC_PHYDYAS<Threshold_FBMC_PHYDYAS_Upper) & (convPSD_FBMC_PHYDYAS>Threshold_FBMC_PHYDYAS_Lower))*(f_FBMC(2)-f_FBMC(1));

    convPSD_FOFDM1=conv(PSD_FOFDM1,ones(size(PSD_FOFDM1)));
    convPSD_FOFDM1(ceil(end/2)+1:end) = [];
    Threshold_FOFDM1_Upper            = (convPSD_FOFDM1(end)-convPSD_FOFDM1(end)*(1-EnergyThreshold)/2);
    Threshold_FOFDM1_Lower            = (convPSD_FOFDM1(end)*(1-EnergyThreshold)/2);
    RequiredBandwidth_FOFDM1          = sum((convPSD_FOFDM1<Threshold_FOFDM1_Upper) & (convPSD_FOFDM1>Threshold_FOFDM1_Lower))*(f_FOFDM1(2)-f_FOFDM1(1));

    convPSD_FOFDM2=conv(PSD_FOFDM2,ones(size(PSD_FOFDM2)));
    convPSD_FOFDM2(ceil(end/2)+1:end) = [];
    Threshold_FOFDM2_Upper            = (convPSD_FOFDM2(end)-convPSD_FOFDM2(end)*(1-EnergyThreshold)/2);
    Threshold_FOFDM2_Lower            = (convPSD_FOFDM2(end)*(1-EnergyThreshold)/2);
    RequiredBandwidth_FOFDM2          = sum((convPSD_FOFDM2<Threshold_FOFDM2_Upper) & (convPSD_FOFDM2>Threshold_FOFDM2_Lower))*(f_FOFDM2(2)-f_FOFDM2(1));

    convPSD_FOFDM3=conv(PSD_FOFDM3,ones(size(PSD_FOFDM3)));
    convPSD_FOFDM3(ceil(end/2)+1:end) = [];
    Threshold_FOFDM3_Upper            = (convPSD_FOFDM3(end)-convPSD_FOFDM3(end)*(1-EnergyThreshold)/2);
    Threshold_FOFDM3_Lower            = (convPSD_FOFDM3(end)*(1-EnergyThreshold)/2);
    RequiredBandwidth_FOFDM3          = sum((convPSD_FOFDM3<Threshold_FOFDM3_Upper) & (convPSD_FOFDM3>Threshold_FOFDM3_Lower))*(f_FOFDM3(2)-f_FOFDM3(1));


    % Calculate the guard frequency F_G
    GuardFrequency_FBMC_Hermite(i_L) = (RequiredBandwidth_FBMC_Hermite - FBMC_Hermite.PHY.SubcarrierSpacing * FBMC_Hermite.Nr.Subcarriers);
    GuardFrequency_FBMC_PHYDYAS(i_L) = (RequiredBandwidth_FBMC_PHYDYAS - FBMC_PHYDYAS.PHY.SubcarrierSpacing * FBMC_PHYDYAS.Nr.Subcarriers);
    GuardFrequency_FOFDM1(i_L)       = (RequiredBandwidth_FOFDM1 - FOFDM1.PHY.SubcarrierSpacing * FOFDM1.Nr.Subcarriers);
    GuardFrequency_FOFDM2(i_L)       = (RequiredBandwidth_FOFDM2 - FOFDM2.PHY.SubcarrierSpacing * FOFDM2.Nr.Subcarriers);
    GuardFrequency_FOFDM3(i_L)       = (RequiredBandwidth_FOFDM3 - FOFDM3.PHY.SubcarrierSpacing * FOFDM3.Nr.Subcarriers);


    % Calculate the time-frequency efficiency
    TimeFrequencyEfficiency_FMBC_Hermite(i_L) = K*L/(RequiredTime_FBMC_Hermite*RequiredBandwidth_FBMC_Hermite); % K for complex symbols! For real symbols K we need to include the factor 0.5
    TimeFrequencyEfficiency_FMBC_PHYDYAS(i_L) = K*L/(RequiredTime_FBMC_PHYDYAS*RequiredBandwidth_FBMC_PHYDYAS); % K for complex symbols! For real symbols K we need to include the factor 0.5
    TimeFrequencyEfficiency_FOFDM1(i_L)       = K*L/(RequiredTime_FOFDM1*RequiredBandwidth_FOFDM1);
    TimeFrequencyEfficiency_FOFDM2(i_L)       = K*L/(RequiredTime_FOFDM2*RequiredBandwidth_FOFDM2);
    TimeFrequencyEfficiency_FOFDM3(i_L)       = K*L/(RequiredTime_FOFDM3*RequiredBandwidth_FOFDM3);

    TimeFrequencyEfficiency_Kinfinity_FMBC_Hermite(i_L) = 0.5*L/(FBMC_Hermite.PHY.TimeSpacing*RequiredBandwidth_FBMC_Hermite);
    TimeFrequencyEfficiency_Kinfinity_FMBC_PHYDYAS(i_L) = 0.5*L/(FBMC_PHYDYAS.PHY.TimeSpacing*RequiredBandwidth_FBMC_PHYDYAS);
    TimeFrequencyEfficiency_Kinfinity_FOFDM1(i_L)       = L/(FOFDM1.PHY.TimeSpacing*RequiredBandwidth_FOFDM1);
    TimeFrequencyEfficiency_Kinfinity_FOFDM2(i_L)       = L/(FOFDM2.PHY.TimeSpacing*RequiredBandwidth_FOFDM2);
    TimeFrequencyEfficiency_Kinfinity_FOFDM3(i_L)       = L/(FOFDM3.PHY.TimeSpacing*RequiredBandwidth_FOFDM3);

    disp([int2str(i_L/length(M_L)*100) '%']);
end

ColorFBMC           = [0 0 1]*0.5;
ColorFBMC_Hermite   = [0 0 1];
ColorFOFDM          = [1 1 0]*0.7;

figure();
plot(M_L,TimeFrequencyEfficiency_FMBC_Hermite,'Color',ColorFBMC_Hermite);hold on;
plot(M_L,TimeFrequencyEfficiency_FMBC_PHYDYAS,'Color',ColorFBMC);hold on;
plot(M_L,TimeFrequencyEfficiency_Kinfinity_FMBC_Hermite,'Color',ColorFBMC_Hermite);hold on;
plot(M_L,TimeFrequencyEfficiency_Kinfinity_FMBC_PHYDYAS,'Color',ColorFBMC);hold on;
% plot(M_L,TimeFrequencyEfficiency_FOFDM1,':','Color',ColorFOFDM);hold on;
% plot(M_L,TimeFrequencyEfficiency_FOFDM2,':','Color',ColorFOFDM);hold on;
% plot(M_L,TimeFrequencyEfficiency_FOFDM3,':','Color',ColorFOFDM);hold on;
plot(M_L,TimeFrequencyEfficiency_Kinfinity_FOFDM1,'Color',ColorFOFDM);hold on;
plot(M_L,TimeFrequencyEfficiency_Kinfinity_FOFDM2,'Color',ColorFOFDM);hold on;
plot(M_L,TimeFrequencyEfficiency_Kinfinity_FOFDM3,'Color',ColorFOFDM);hold on;
xlim([0 max(M_L)]);
xlabel('Number of Subcarriers L');
ylabel('Time-Frequency Efficiency $\rho$');
title('See Fig. 9');


