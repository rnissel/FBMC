% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script plots the power spectal density of different modulation
% techniques:
%       1) OFDM with CP (worst spectral behaviour)
%       2) FBMC (best spectral properties, complex orthogonality is replaced by real orthogonality)
%       3) WOLA (windowed OFDM. The windowing is done at TX and RX)
%       4) FOFDM (filtered OFDM, sinc + Hann, filtering at TX and RX)
%       5) UFMC (filtered OFDM, subband-wise filtering, Dolph-Chebyshev window, cyclic prefix and zero padding are supported, filtering at TX and RX)
% See Figure 2 in the paper

clear; close all;

L                           = 24;                                       % Number of subcarriers
K_OFDM                      = 10;                                       % Number of OFDM symbols in time
K_FBMC                      = 105;                                      % For FBMC we use a higher number of time-symbols so that the frequency resolution is better. Without that, the figure would look ugly.

NrRepetitions               = 200;                                      % Number repetitions for the simulation
QAM_ModulationOrder         = 16;                                       % Modulation order, 4,16,64,...

Simulate                    = true;                                     % Perform also a simulation in order to check the theoretical Power Spectral Density (PSD)

SubcarrierSpacing           = 15e3;                                     % Subcarrier spacing (15kHz, same as LTE)
SamplingRate                = SubcarrierSpacing*L*14;                   % We need oversampling (14) in order to see out of band emissions. Furthermore 14 fits the CP length of OFDM.

OverlappingFactor           = 4;                                        % Overlapping factor for FBMC. For an overlapping factor of 4 and the PHYDYAS pulse, we see a rectangular filering effect     
IntermediateSubcarrier      = 50;                                       % Shift frequency (for presentation purposes)
UFMC_ZeroPadding            = false;                                    % If true, apply zero padding instead of a cyclic prefix. However, zero padding leads to zero crossing which then dominate the whole picture => for presentation purposes we also use CP in UFMC  

% Parameterset 1
TF_FilteredAndWindowedOFDM  = 1.09;                                     % Time frequency spacing for UFMC and FOFDM. Chosen so that the interference from filtering is small (the SIR is approximately 65dB).
FilterLengthTXandRX_UFMC    = 1/(14*SubcarrierSpacing);                 % 4.7µs, same as the CP in LTE. However, we need filtering at TX and RX => Total length of the filtering is two times larger
FilterLengthTXandRX_FOFDM   = 1.18*1/(14*SubcarrierSpacing);            % Increased by 1.2 because FOFDM is more robust than UFMC  to a larger filter length

% Parameterset 2
% TF_FilteredAndWindowedOFDM = 1.28;                                    % Time frequency spacing for UFMC and FOFDM. Chosen so that the interference from filtering is small (the SIR is approximately 65dB).
% FilterLengthTXandRX_UFMC = 3/(14*SubcarrierSpacing);                  % 3 times as high as the CP in LTE to get a smaller OOB. 
% FilterLengthTXandRX_FOFDM = 1.33*3/(14*SubcarrierSpacing);            % Increased by 1.35 because FOFDM is more robust than UFMC to a larger filter length

% Parameterset 3 (perfect orthogonal)
% TF_FilteredAndWindowedOFDM = 1+2/(14);
% FilterLengthTXandRX_UFMC = 1/(14*SubcarrierSpacing);        
% FilterLengthTXandRX_FOFDM = 1/(14*SubcarrierSpacing);

% Parameterset 4 (LTE like, interference ~ 55)
% TF_FilteredAndWindowedOFDM  = 1+1/(14);
% FilterLengthTXandRX_UFMC    = 1/(14*SubcarrierSpacing);        
% FilterLengthTXandRX_FOFDM   = 1.1*1/(14*SubcarrierSpacing);

ColorFBMC   = [0 0 1]*0.5;
ColorOFDM   = [1 0 0];
ColorWOLA   = [1 0 1];
ColorFOFDM  = [1 1 0]*0.7;
ColorUFMC   = [0 1 0]*0.5;

CP_LengthFilter = (TF_FilteredAndWindowedOFDM-1)*1/(SubcarrierSpacing);

%% FBMC Object
FBMC = Modulation.FBMC(...
    L,...                                                       % Number subcarriers
    K_FBMC,...                                                  % Number FBMC symbols (determines the frequency resolution!)
    SubcarrierSpacing,...                                       % Subcarrier spacing (Hz)
    SamplingRate,...                                            % Sampling rate (Samples/s)
    SubcarrierSpacing*IntermediateSubcarrier,...                % Intermediate frequency first subcarrier (Hz)
    false,...                                                   % Transmit real valued signal
    'PHYDYAS-OQAM',...                                          % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    OverlappingFactor, ...                                      % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                      % Initial phase shift
    true ...                                                    % Polyphase implementation
    );
%% OFDM Object
OFDM = Modulation.OFDM(...
    L,...                                                       % Number subcarriers
    K_OFDM,...                                                  % Number OFDM Symbols
    SubcarrierSpacing,...                                       % Subcarrier spacing (Hz)
    SamplingRate,...                                            % Sampling rate (Samples/s)
    SubcarrierSpacing*IntermediateSubcarrier,...                % Intermediate frequency first subcarrier (Hz)
    false,...                                                   % Transmit real valued signal
    1/(14*SubcarrierSpacing), ...                               % Cyclic prefix length (s) 
    0 ...                                                       % Zero guard length (s)
    );

%% Windowed OFDM (WOLA)
WOLA = Modulation.WOLA(...
    L,...                                                       % Number subcarriers
    K_OFDM,...                                                  % Number OFDM Symbols
    SubcarrierSpacing,...                                       % Subcarrier spacing (Hz)
    SamplingRate,...                                            % Sampling rate (Samples/s)
    SubcarrierSpacing*IntermediateSubcarrier,...                % Intermediate frequency first subcarrier (Hz)
    false,...                                                   % Transmit real valued signal
    0, ...                                                      % Cyclic prefix length (s) 
    0, ...                                                      % Zero guard length (s)
    CP_LengthFilter/2, ...                                      % Window overlap length at TX (s)
    CP_LengthFilter/2 ...                                       % Window overlap length at RX (s)
    );

%% Filtered OFDM 
FOFDM = Modulation.FOFDM(...
    L,...                                                       % Number subcarriers
    K_OFDM,...                                                  % Number OFDM Symbols
    SubcarrierSpacing,...                                       % Subcarrier spacing (Hz)
    SamplingRate,...                                            % Sampling rate (Samples/s)
    SubcarrierSpacing*IntermediateSubcarrier,...                % Intermediate frequency first subcarrier (Hz)
    false,...                                                   % Transmit real valued signal
    0, ...                                                      % Cyclic prefix length (s) 
    0, ...                                                      % Zero guard length (s)
    FilterLengthTXandRX_FOFDM, ...                              % Filter length at TX (s)
    FilterLengthTXandRX_FOFDM, ...                              % Filter length at RX (s)
    CP_LengthFilter ...                                         % Addional cyclic prefix for the filtering (s)
);

%% UFMC
UFMC = Modulation.UFMC(...
    L,...                                                       % Number subcarriers
    K_OFDM,...                                                  % Number OFDM Symbols
    SubcarrierSpacing,...                                       % Subcarrier spacing (Hz)
    SamplingRate,...                                            % Sampling rate (Samples/s)
    SubcarrierSpacing*IntermediateSubcarrier,...                % Intermediate frequency first subcarrier (Hz)
    false,...                                                   % Transmit real valued signal
    0, ...                                                      % Cyclic prefix length (s) 
    0, ...                                                      % Zero guard length (s)
    FilterLengthTXandRX_UFMC, ...                               % Filter length at TX (s)
    FilterLengthTXandRX_UFMC, ...                               % Filter length at RX (s)
    CP_LengthFilter, ...                                        % Addional cyclic prefix for the filtering (s)
    UFMC_ZeroPadding ...                                        % TRUE for zero padding; FALSE for a conventional cyclic prefix
);


%% Get TX and RX matrices to calculate SIR for FOFDM and UFMC due to ISI/ICI
G_TX_FOFDM = sparse(FOFDM.GetTXMatrix);
G_RX_FOFDM = sparse(FOFDM.GetRXMatrix);

G_TX_UFMC = sparse(UFMC.GetTXMatrix);
G_RX_UFMC = sparse(UFMC.GetRXMatrix);

D_FOFDM = G_RX_FOFDM*G_TX_FOFDM;
D_UFMC  = G_RX_UFMC*G_TX_UFMC;

% Inband Signal to interference ratio due to filtering (the filter length is longer than the CP to improve filter characteristics)
SIR_dB_FOFDM = 10*log10(sum(abs(diag(D_FOFDM)).^2)/(sum(sum(abs(D_FOFDM-diag(diag(D_FOFDM))).^2,2),1)));
SIR_dB_UFMC  = 10*log10(sum(abs(diag(D_UFMC)).^2)/(sum(sum(abs(D_UFMC-diag(diag(D_UFMC))).^2,2),1)));


%% Information
fprintf('===================================================\n');
fprintf('               |(complex)TF-Spacing| Bandwidth(LF)| \n');
fprintf('OFDM (with CP) |%17.2f  |%8.2f MHz  | \n', OFDM.PHY.TimeSpacing*OFDM.PHY.SubcarrierSpacing , OFDM.PHY.SubcarrierSpacing*OFDM.Nr.Subcarriers/1e6);
fprintf('FBMC           |%17.2f  |%8.2f MHz  | \n', FBMC.PHY.TimeSpacing*FBMC.PHY.SubcarrierSpacing*2 , FBMC.PHY.SubcarrierSpacing*FBMC.Nr.Subcarriers/1e6);
fprintf('WOLA           |%17.2f  |%8.2f MHz  | \n', WOLA.PHY.TimeSpacing*WOLA.PHY.SubcarrierSpacing , WOLA.PHY.SubcarrierSpacing*WOLA.Nr.Subcarriers/1e6);
fprintf('FOFDM          |%17.2f  |%8.2f MHz  | \n', FOFDM.PHY.TimeSpacing*FOFDM.PHY.SubcarrierSpacing , FOFDM.PHY.SubcarrierSpacing*FOFDM.Nr.Subcarriers/1e6);
fprintf('UFMC           |%17.2f  |%8.2f MHz  | \n', UFMC.PHY.TimeSpacing*UFMC.PHY.SubcarrierSpacing , UFMC.PHY.SubcarrierSpacing*UFMC.Nr.Subcarriers/1e6);
fprintf('===================================================\n');
fprintf('SIR [dB] for FOFDM: %2.2f   \n', full(SIR_dB_FOFDM));
fprintf('SIR [dB] for UFMC:  %2.2f   \n', full(SIR_dB_UFMC));
fprintf('SIR [dB] for FBMC:  %2.2f   \n', FBMC.GetSIRdBDoublyFlat);  % Interference due to imperfect prototype filter


%% Calculate the Power Spectral Density (PSD)
[PSD_FBMC_Theory,f_FBMC]   = FBMC.PlotPowerSpectralDensityUncorrelatedData;
[PSD_OFDM_Theory,f_OFDM]   = OFDM.PlotPowerSpectralDensityUncorrelatedData;
[PSD_WOLA_Theory,f_WOLA]   = WOLA.PlotPowerSpectralDensityUncorrelatedData;
[PSD_FOFDM_Theory,f_FOFDM] = FOFDM.PlotPowerSpectralDensityUncorrelatedData;
[PSD_UFMC_Theory,f_UFMC]   = UFMC.PlotPowerSpectralDensityUncorrelatedData;

% Normalize over energy
PSD_FBMC_Theory  = PSD_FBMC_Theory/(sum(PSD_FBMC_Theory)*SamplingRate/FBMC.Nr.SamplesTotal);
PSD_OFDM_Theory  = PSD_OFDM_Theory/(sum(PSD_OFDM_Theory)*SamplingRate/OFDM.Nr.SamplesTotal);
PSD_WOLA_Theory  = PSD_WOLA_Theory/(sum(PSD_WOLA_Theory)*SamplingRate/WOLA.Nr.SamplesTotal);
PSD_FOFDM_Theory = PSD_FOFDM_Theory/(sum(PSD_FOFDM_Theory)*SamplingRate/FOFDM.Nr.SamplesTotal);
PSD_UFMC_Theory  = PSD_UFMC_Theory/(sum(PSD_UFMC_Theory)*SamplingRate/UFMC.Nr.SamplesTotal);

% Normalize to 0dB
NormalizationFactor_Theory = max([PSD_FBMC_Theory]);
PSD_FBMC_Theory  = PSD_FBMC_Theory/NormalizationFactor_Theory;
PSD_OFDM_Theory  = PSD_OFDM_Theory/NormalizationFactor_Theory;
PSD_WOLA_Theory  = PSD_WOLA_Theory/NormalizationFactor_Theory;
PSD_FOFDM_Theory = PSD_FOFDM_Theory/NormalizationFactor_Theory;
PSD_UFMC_Theory  = PSD_UFMC_Theory/NormalizationFactor_Theory;


figure();
plot([-12 -12], [-150,20],'color', [0.8 0.8 0.8]);
hold on;
plot([12 12], [-150,20],'color', [0.8 0.8 0.8]);
plot(f_OFDM/OFDM.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_OFDM_Theory]),'Color', ColorOFDM);
plot(f_WOLA/WOLA.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_WOLA_Theory]),'Color', ColorWOLA);
plot(f_UFMC/UFMC.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_UFMC_Theory]),'Color',ColorUFMC);
plot(f_FOFDM/FOFDM.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_FOFDM_Theory]),'Color',ColorFOFDM);
plot(f_FBMC/FBMC.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_FBMC_Theory]),'Color', ColorFBMC);
ylim([-100 4]);
xlim([-50 50]);
xlabel('Normalized Frequency, f/F');
ylabel('Power Spectral Density [dB]');
set(gca,'XTick',[-50:10:50]);


if Simulate
%% Simulate to Check Theory
tic
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
PAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
PSD_FBMC_Simulation  = zeros(FBMC.Nr.SamplesTotal,1);
PSD_OFDM_Simulation  = zeros(OFDM.Nr.SamplesTotal,1);
PSD_WOLA_Simulation  = zeros(WOLA.Nr.SamplesTotal,1);
PSD_FOFDM_Simulation = zeros(FOFDM.Nr.SamplesTotal,1);
PSD_UFMC_Simulation  = zeros(UFMC.Nr.SamplesTotal,1);
for i_rep = 1:NrRepetitions
    x_OFDM = QAM.SymbolMapping(randi(QAM_ModulationOrder,L,K_OFDM));
    x_PAM  = PAM.SymbolMapping(randi(sqrt(QAM_ModulationOrder),L,K_FBMC));

    s_FBMC  = FBMC.Modulation(x_PAM);
    s_OFDM  = OFDM.Modulation(x_OFDM);
    s_WOLA  = WOLA.Modulation(x_OFDM);
    s_FOFDM = FOFDM.Modulation(x_OFDM);
    s_UFMC  = UFMC.Modulation(x_OFDM);
    
    PSD_FBMC_Simulation  = PSD_FBMC_Simulation + abs(fft(s_FBMC)).^2;   
    PSD_OFDM_Simulation  = PSD_OFDM_Simulation + abs(fft(s_OFDM)).^2;
    PSD_WOLA_Simulation  = PSD_WOLA_Simulation + abs(fft(s_WOLA)).^2;
    PSD_FOFDM_Simulation = PSD_FOFDM_Simulation + abs(fft(s_FOFDM)).^2;
    PSD_UFMC_Simulation  = PSD_UFMC_Simulation + abs(fft(s_UFMC)).^2;
    
    
    TimePassed = toc;
    if mod(i_rep,100)==0
        disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes']);
    end
end
% Normalize over energy
PSD_FBMC_Simulation = PSD_FBMC_Simulation/(sum(PSD_FBMC_Simulation)*SamplingRate/FBMC.Nr.SamplesTotal);
PSD_OFDM_Simulation = PSD_OFDM_Simulation/(sum(PSD_OFDM_Simulation)*SamplingRate/OFDM.Nr.SamplesTotal);
PSD_WOLA_Simulation = PSD_WOLA_Simulation/(sum(PSD_WOLA_Simulation)*SamplingRate/WOLA.Nr.SamplesTotal);
PSD_FOFDM_Simulation = PSD_FOFDM_Simulation/(sum(PSD_FOFDM_Simulation)*SamplingRate/FOFDM.Nr.SamplesTotal);
PSD_UFMC_Simulation = PSD_UFMC_Simulation/(sum(PSD_UFMC_Simulation)*SamplingRate/UFMC.Nr.SamplesTotal);

% LS Estimation to normalize to 0dB
NormalizationFactor_Simulation = 1./mean((PSD_FBMC_Simulation'*PSD_FBMC_Theory)./(PSD_FBMC_Simulation'*PSD_FBMC_Simulation));
PSD_FBMC_Simulation = PSD_FBMC_Simulation/NormalizationFactor_Simulation;
PSD_OFDM_Simulation = PSD_OFDM_Simulation/NormalizationFactor_Simulation;
PSD_WOLA_Simulation = PSD_WOLA_Simulation/NormalizationFactor_Simulation;
PSD_FOFDM_Simulation = PSD_FOFDM_Simulation/NormalizationFactor_Simulation;
PSD_UFMC_Simulation = PSD_UFMC_Simulation/NormalizationFactor_Simulation;

figure();
plot(f_OFDM/OFDM.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_OFDM_Theory PSD_OFDM_Simulation]),'Color', ColorOFDM);
hold on;
plot(f_WOLA/WOLA.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_WOLA_Theory PSD_WOLA_Simulation]),'Color', ColorWOLA);
plot(f_UFMC/UFMC.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_UFMC_Theory PSD_UFMC_Simulation]),'Color', ColorUFMC);
plot(f_FOFDM/FOFDM.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_FOFDM_Theory PSD_FOFDM_Simulation]),'Color', ColorFOFDM);
plot(f_FBMC/FBMC.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([PSD_FBMC_Theory PSD_FBMC_Simulation]),'Color', ColorFBMC);
plot([-12 -12], [-150,20],'color', [0.8 0.8 0.8]);
plot([12 12], [-150,20],'color', [0.8 0.8 0.8]);
ylim([-100 4]);
xlim([-50 50]);
xlabel('Normalized Frequency, f/F');
ylabel('Power Spectral Density [dB]');
set(gca,'XTick',[-50:10:50]);
title('Theory vs Simulation');

end

