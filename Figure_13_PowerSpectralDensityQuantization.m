% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script plots the power spectal density of FBMC and OFDM in case of
% quantization and clipping.
% See Figure 13 in the paper

clear; close all;

% Important parameters for the quantization
M_BitResolution = [4 8 12 16];                              % DAC resolution: 4-bit, 8-bit, 12-bit, 16-bit
M_Clipping_dB   = [8 9 11 12];                              % The signal has a mean power of one. Values larger(smaller) than Clipping_dB are set to the maximum (minimum) value and affect the Power Spectral Density (PSD)


% additional parameters
NrRepetitions_PSD   = 100; %10000                           % Number repetitions for the simulation
K_FBMC              = 105;                                  % Number of FBMC symbols
K_OFDM              = 10;                                   % Number of OFDM symbols
L                   = 24;                                   % Number of subcarriers
QAM_ModulationOrder = 64;                                   % Modulation order, 4,16,64,...
SubcarrierSpacing   = 15e3;                                 % Subcarrier spacing (15kHz, same as LTE)
SamplingRate        = SubcarrierSpacing*1024;               % A high oversampling factor reduces the quantization effect. We therfore chose 24*42.66=1024 instead of 24*14.

OverlappingFactor = 4;                                      % Overlapping factor for FBMC. For an overlapping factor of 4 and the PHYDYAS pulse, we see a rectangular filering effect
IntermediateSubcarrier = 50;                                % Shift frequency (for presentation purposes)

ColorFBMC   = [0 0 1]*0.5;
ColorOFDM   = [1 0 0];

%% FBMC Object
FBMC = Modulation.FBMC(...
    L,...                                                   % Number subcarriers
    K_FBMC,...                                              % Number FBMC symbols
    SubcarrierSpacing,...                                   % Subcarrier spacing (Hz)
    SamplingRate,...                                        % Sampling rate (Samples/s)
    SubcarrierSpacing*IntermediateSubcarrier,...            % Intermediate frequency first subcarrier (Hz)
    false,...                                               % Transmit real valued signal
    'PHYDYAS-OQAM',...                                      % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    OverlappingFactor, ...                                  % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                  % Initial phase shift
    true ...                                                % Polyphase implementation
    );
%% OFDM Object
OFDM = Modulation.OFDM(...
    L,...                                                   % Number subcarriers
    K_OFDM,...                                              % Number OFDM Symbols
    SubcarrierSpacing,...                                   % Subcarrier spacing (Hz)
    SamplingRate,...                                        % Sampling rate (Samples/s)
    SubcarrierSpacing*IntermediateSubcarrier,...            % Intermediate frequency first subcarrier (Hz)
    false,...                                               % Transmit real valued signal
    1/(14*SubcarrierSpacing), ...                           % Cyclic prefix length (s) 
    (OverlappingFactor-1/2)*1/SubcarrierSpacing*1/2 ...     % Zero guard length (s)
    );

%% Modulation Object
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
PAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');

%% Power Spectral Density
for i_BitResolution = 1:size(M_BitResolution,2)
BitResolution   = M_BitResolution(i_BitResolution);
Clipping_dB     = M_Clipping_dB(i_BitResolution);
disp('===========================================');
disp(['Bit resolution: ' int2str(BitResolution)]);
disp(['Clipping: ' int2str(Clipping_dB) 'dB']);
% PSD...Power Spectral Density
PSD_OFDM_Simulation = zeros(OFDM.Nr.SamplesTotal,1);
PSD_FBMC_Simulation = zeros(FBMC.Nr.SamplesTotal,1);
Scaling = 2^(BitResolution-1)/10^(Clipping_dB/20);
tic
for i_rep = 1:NrRepetitions_PSD
    x_OFDM  = QAM.SymbolMapping(randi(QAM_ModulationOrder,OFDM.Nr.Subcarriers,OFDM.Nr.MCSymbols));
    x_PAM   = PAM.SymbolMapping(randi(sqrt(QAM_ModulationOrder),FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols));   
    
    s_OFDM  = OFDM.Modulation(x_OFDM);
    s_FBMC  = FBMC.Modulation(x_PAM);
    
    % Scale and quantization
    s_OFDM  = round(s_OFDM*Scaling);
    s_FBMC  = round(s_FBMC*Scaling);
    
    % Clipping OFDM
    Real_s_OFDM                                         = real(s_OFDM);
    Imag_s_OFDM                                         = imag(s_OFDM);    
    NumberOfClippedSamples_OFDM(i_rep)                  = sum(Real_s_OFDM>(2^(BitResolution-1)-1))+sum(Real_s_OFDM<-(2^(BitResolution-1))); % just for information
    Real_s_OFDM(Real_s_OFDM > (2^(BitResolution-1)-1))  = (2^(BitResolution-1)-1);
    Real_s_OFDM(Real_s_OFDM <-(2^(BitResolution-1)))    = -(2^(BitResolution-1));    
    Imag_s_OFDM(Imag_s_OFDM > (2^(BitResolution-1)-1))  = (2^(BitResolution-1)-1);
    Imag_s_OFDM(Imag_s_OFDM <-(2^(BitResolution-1)))    = -(2^(BitResolution-1));       
    s_OFDM                                              =  Real_s_OFDM+1j*Imag_s_OFDM;
    
    % Clipping FBMC
    Real_s_FBMC                                         = real(s_FBMC); 
    Imag_s_FBMC                                         = imag(s_FBMC);   
    NumberOfClippedSamples_FBMC(i_rep)                  = sum(Real_s_FBMC>(2^(BitResolution-1)-1))+sum(Real_s_FBMC<-(2^(BitResolution-1)));
    Real_s_FBMC(Real_s_FBMC > (2^(BitResolution-1)-1))  = (2^(BitResolution-1)-1);
    Real_s_FBMC(Real_s_FBMC <-(2^(BitResolution-1)))    = -(2^(BitResolution-1));    
    Imag_s_FBMC(Imag_s_FBMC > (2^(BitResolution-1)-1))	= (2^(BitResolution-1)-1);
    Imag_s_FBMC(Imag_s_FBMC <-(2^(BitResolution-1)))    = -(2^(BitResolution-1));       
    s_FBMC                                              =  Real_s_FBMC+1j*Imag_s_FBMC;    
    
    % Calculate PSD
    PSD_OFDM_Simulation = PSD_OFDM_Simulation + abs(fft(s_OFDM)).^2;
    PSD_FBMC_Simulation = PSD_FBMC_Simulation + abs(fft(s_FBMC)).^2;
    
    % Calculate Interference Power
    PI_OFDM(i_rep) = mean(mean(abs(OFDM.Demodulation(s_OFDM)/Scaling-x_OFDM).^2));
    PI_FBMC(i_rep) = mean(mean(abs(real(FBMC.Demodulation(s_FBMC))/Scaling-x_PAM).^2));
        
    TimePassed = toc;
    if mod(i_rep,100)==0
        disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions_PSD) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions_PSD-i_rep)/60) 'minutes']);
    end
end

% Normalize over energy
PSD_FBMC_Simulation = PSD_FBMC_Simulation/(sum(PSD_FBMC_Simulation)*SamplingRate/FBMC.Nr.SamplesTotal);
PSD_OFDM_Simulation = PSD_OFDM_Simulation/(sum(PSD_OFDM_Simulation)*SamplingRate/OFDM.Nr.SamplesTotal);

NormFactor = max([PSD_FBMC_Simulation(:)]);
PSD_OFDM_Simulation = PSD_OFDM_Simulation/NormFactor;
PSD_FBMC_Simulation = PSD_FBMC_Simulation/NormFactor;

disp('Number of clipped samples: OFDM vs FBMC');
disp([sum(NumberOfClippedSamples_OFDM)  sum(NumberOfClippedSamples_FBMC)]);

M_SIR_OFDM(i_BitResolution) = -10*log10(mean(PI_OFDM));
M_SIR_FBMC(i_BitResolution) = -10*log10(mean(PI_FBMC));

M_PSD_OFDM_Simulation(:,i_BitResolution) = PSD_OFDM_Simulation;
M_PSD_FBMC_Simulation(:,i_BitResolution) = PSD_FBMC_Simulation;

end


f_OFDM = (0:(size(PSD_OFDM_Simulation,1)-1))*1/(size(PSD_OFDM_Simulation,1)*OFDM.PHY.dt);
f_FBMC = (0:(size(PSD_FBMC_Simulation,1)-1))*1/(size(PSD_FBMC_Simulation,1)*FBMC.PHY.dt);

f_Normalized_OFDM = f_OFDM/OFDM.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1;
f_Normalized_FBMC = f_FBMC/FBMC.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1;


figure();
plot([-12 -12], [-150,20],'color', [0.8 0.8 0.8]);
hold on;
plot([12 12], [-150,20],'color', [0.8 0.8 0.8]);
plot(f_OFDM/OFDM.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([M_PSD_OFDM_Simulation(:,end).']),'Color', ColorOFDM);
plot(f_FBMC/FBMC.PHY.SubcarrierSpacing-IntermediateSubcarrier-ceil(L/2)+1/2,10*log10([M_PSD_FBMC_Simulation.']),'Color', ColorFBMC);
ylim([-100 4]);
xlim([-50 50]);
xlabel('Normalized Frequency, $f/F$');
ylabel('Power Spectral Density [dB]');
set(gca,'XTick',[-50:10:50]);

disp('SNR [dB] for OFDM;');
disp(M_SIR_OFDM);
disp('SNR [dB] for FBMC;'); % Increasing the overlapping factor increases the SIR bound
disp(M_SIR_FBMC);





