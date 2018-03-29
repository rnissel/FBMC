% =====================================================================    
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =====================================================================    
% This script plots the Signal-to-Interference Ratio (SIR) for an optimal
% subcarrier spacing. The SIR values, depending on velocitiy and subcarrier
% spacing, are precalulated in the folder "./Results". This script then
% combines those pre-calulated SIR values.

clear; close all;

PlotAddionalResults = false;                            % If true, plot additional information.    



%% Plot Pedestrian A  
ColorFBMC   = [0 0 1]*0.5; % Color for PHYDYAS pulse
M_Velocity_kmh = [0:2:8 10:10:500];        
M_SubcarrierSpacingInterpolated = 1e3*[5:0.1:180];    
NrMethods = 4;
PowerDelayProfile = 'PedestrianA';
Method{1} = 'OFDM';
PrototypeFilter{1} = '';        
Method{2} = 'FBMC-OQAM';
PrototypeFilter{2} = 'Hermite';
Method{3} = 'FBMC-OQAM';
PrototypeFilter{3} = 'PHYDYAS';  
Method{4} = 'CP-OFDM-007LTE';
PrototypeFilter{4} = '';           
  
% Load results and find optimum (using interpolation to obtain smoother curves)
for i_method = 1:NrMethods
    for i_velocity = 1:length(M_Velocity_kmh)   
        Velocity_kmh = M_Velocity_kmh(i_velocity);
        FileName = ['.\Results\' Method{i_method} '_' PrototypeFilter{i_method} '_' PowerDelayProfile int2str(Velocity_kmh) '.mat'];
        L1 = load(FileName);

        % Needed for interpolation: remove double subcarriers
        IndexSameSubcarrierSpacing = find(diff(L1.M_RealSubcarrierSpacing)==0);
        L1.M_RealSubcarrierSpacing(IndexSameSubcarrierSpacing) = [];
        L1.SIR_Theory(IndexSameSubcarrierSpacing) = [];

        % Interpolate to get smoother curves
        M_SIR_Interpolated(:,i_velocity) = interp1(L1.M_RealSubcarrierSpacing,L1.SIR_Theory,M_SubcarrierSpacingInterpolated,'spline',nan);

        if strcmp(Method{i_method},'CP-OFDM-007LTE') % To check if theoretical equations are correct (interpolation is wrong in case of CP-OFDM!)
            M_SIR007_Interpolated_True(:,i_velocity) = L1.SIR_Theory;
        end
    end
    M_SIR_Interpolated_ALL(:,:,i_method) = M_SIR_Interpolated;

    [SIRoptimal,IndexSIRoptimal] = max(M_SIR_Interpolated);
    SubcarrierSpacingOptimal = M_SubcarrierSpacingInterpolated(IndexSIRoptimal);

    M_SIR_15kHZ(:,i_method) = M_SIR_Interpolated(M_SubcarrierSpacingInterpolated==15e3,:);
    M_SIRoptimal(:,i_method) = SIRoptimal.';
    M_SubcarrierSpacingOptimal(:,i_method) = SubcarrierSpacingOptimal.';
end
M_SIR_15kHZ(:,4)=nan;               % Values are wrong for CP-OFDM due to discontinuity
M_SIRoptimal(:,4)=nan;
M_SubcarrierSpacingOptimal(:,4)=nan;


% Calculate theoretical SIR for OFDM (no ISI)
Kappa           = 1/14;             % LTE efficiency loss
TCP             = 410e-9;           % maximum delay of Pedestrian A
F               = Kappa/TCP;        % optimal subcarrier spacing (no ISI)
M_DopplerShift  = M_Velocity_kmh/3.6*2.5e9/2.998e8;
P_S_CPOFDM      = hypergeom([1/2],[3/2 2],-(M_DopplerShift/F*pi).^2);
P_I_CPOFDM      = 1-P_S_CPOFDM;
SIR_CPOFDM      = 10*log10(P_S_CPOFDM./P_I_CPOFDM);


figure(11);
plot(M_Velocity_kmh,M_SIRoptimal(:,1),':','color', [0.5 0.5 0]);
hold on;
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,1)>5e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,1)>5e3,1),'color', [0.5 0.5 0]);
plot(M_Velocity_kmh,M_SIRoptimal(:,2),': blue');        
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,2)>5e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,2)>5e3,2),'blue');        
plot(M_Velocity_kmh,M_SIRoptimal(:,3),'color',ColorFBMC); 
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,3)>5e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,3)>5e3,3),'color',ColorFBMC);  
plot(M_Velocity_kmh,SIR_CPOFDM,'red');

xlabel('Velocity [km/h]');
ylabel('Signal-to-Interference Ratio [dB]');
ylim([0 50]);
p1 = plot(nan,nan,': black');
p2 = plot(nan,nan,'black');
legend([p1 p2],{'Subcarrier Spacing, F=5kHz','Subcarrier Spacing, F>5kHz'},'Location','SouthWest');
 
if PlotAddionalResults
    % Optimal subcarrier spacing according to the rule of thumb   
    RMS_Doppler = 1/sqrt(2)*M_DopplerShift;
    RMS_Tau = 46e-9;

    ApproxOptimalSubcarrierSpacing_PHYDYAS = sqrt(0.2745/0.328*RMS_Doppler/RMS_Tau);
    ApproxOptimalSubcarrierSpacing_Hermite = sqrt(0.2015/0.403*RMS_Doppler/RMS_Tau);

    figure(12);
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,1)>5e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,1)>5e3,1)/1e3,'red');
    hold on;
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,2)>15e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,2)>15e3,2)/1e3,'blue');        
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_Hermite/1e3,': blue');        
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,3)>15e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,3)>15e3,3)/1e3,'color',ColorFBMC);
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_PHYDYAS/1e3,':','color',ColorFBMC);        
    xlabel('Velocity [km/h], f_c=2.5GHz$');
    ylabel('Optimal Subcarrier Spacing [kHz]');
    title('Pedestrian A');
    p1 = plot(nan,nan,'black');
    p2 = plot(nan,nan,': black');
    legend([p1 p2],{'Exhaustive search','Approximation'},'Location','NorthWest');

    % SIR for approximate optimal subcarrier spacing
    rep1 = length(ApproxOptimalSubcarrierSpacing_PHYDYAS);
    rep2 = length(M_SubcarrierSpacingInterpolated);
    [~,index_PHYDYAS] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_PHYDYAS,[rep2,1])),[],1);
    [~,index_Hermite] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_Hermite,[rep2,1])),[],1);

    M_SIR_Interploated_PHYDYAS = M_SIR_Interpolated_ALL(:,:,3);
    M_SIR_Interploated_Hermite = M_SIR_Interpolated_ALL(:,:,2);
    for i_velocity = 1:length(M_Velocity_kmh)
        M_SIR_ApproxOpt_PHYDYAS(i_velocity) = M_SIR_Interploated_PHYDYAS(index_PHYDYAS(i_velocity),i_velocity);
        M_SIR_ApproxOpt_Hermite(i_velocity) = M_SIR_Interploated_Hermite(index_Hermite(i_velocity),i_velocity);

    end

    figure(13);plot(M_SIRoptimal(:,3)-M_SIR_ApproxOpt_PHYDYAS.');
    hold on;
    plot(M_SIRoptimal(:,2)-M_SIR_ApproxOpt_Hermite.');
    xlabel('Velocity [km/h], f_c=2.5GHz');
    ylabel('SIR Error in dB!!! between approx. and exhaustive search');
    title('Pedestrian A');
    % SIR of theory vs simulations
    % plot(M_Velocity_kmh,[max(M_SIR007_Interpolated_True,[],1).' SIR_CPOFDM.']);
end




%% Plot Vehicular A
clear M_SIR_Interpolated_ALL M_SIRoptimal M_SIR_Interpolated M_SIR_15kHZ M_SubcarrierSpacingOptimal M_SIR007_Interpolated_True;
M_Velocity_kmh = [0 2:2:8 10:10:500];        
M_SubcarrierSpacingInterpolated = 1e3*[5:0.1:60];
NrMethods = 4;
PowerDelayProfile = 'VehicularA';
Method{1} = 'OFDM';
PrototypeFilter{1} = '';        
Method{2} = 'FBMC-OQAM';
PrototypeFilter{2} = 'Hermite';
Method{3} = 'FBMC-OQAM';
PrototypeFilter{3} = 'PHYDYAS';
Method{4} = 'CP-OFDM-007LTE';
PrototypeFilter{4} = '';             

% Load results and find optimum (using interpolation to obtain smoother curves)
for i_method = 1:NrMethods
    for i_velocity = 1:length(M_Velocity_kmh)   
        Velocity_kmh = M_Velocity_kmh(i_velocity);
        FileName = ['.\Results\' Method{i_method} '_' PrototypeFilter{i_method} '_' PowerDelayProfile int2str(Velocity_kmh) '.mat'];
        L1 = load(FileName);

        % Needed for interpolation: remove double subcarriers
        IndexSameSubcarrierSpacing = find(diff(L1.M_RealSubcarrierSpacing)==0);
        L1.M_RealSubcarrierSpacing(IndexSameSubcarrierSpacing) = [];
        L1.SIR_Theory(IndexSameSubcarrierSpacing) = [];

        % Interpolate to get smoother curves
        M_SIR_Interpolated(:,i_velocity) = interp1(L1.M_RealSubcarrierSpacing,L1.SIR_Theory,M_SubcarrierSpacingInterpolated,'spline',nan);

        if strcmp(Method{i_method},'CP-OFDM-007LTE') % To check if theoretical equations are correct (interpolation is wrong in case of CP-OFDM!)
            M_SIR007_Interpolated_True(:,i_velocity) = L1.SIR_Theory;
        end
    end
    M_SIR_Interpolated_ALL(:,:,i_method) = M_SIR_Interpolated;


    [SIRoptimal,IndexSIRoptimal] = max(M_SIR_Interpolated);
    SubcarrierSpacingOptimal = M_SubcarrierSpacingInterpolated(IndexSIRoptimal);

    M_SIR_15kHZ(:,i_method) = M_SIR_Interpolated(M_SubcarrierSpacingInterpolated==15e3,:);
    M_SIRoptimal(:,i_method) = SIRoptimal.';
    M_SubcarrierSpacingOptimal(:,i_method) = SubcarrierSpacingOptimal.';
end
M_SIR_15kHZ(:,4)=nan; % Values are wrong for CP-OFDM due to discontinuity
M_SIRoptimal(:,4)=nan;
M_SubcarrierSpacingOptimal(:,4)=nan;


% Calculate Theoretical SIR for OFDM (no ISI)
Kappa = 1/14;   % LTE efficiency loss
TCP = 2510e-9;   % maximum delay of Pedestrian A
F = Kappa/TCP;  %optimal subcarrier spacing (no ISI)
M_DopplerShift = M_Velocity_kmh/3.6*2.5e9/2.998e8;
P_S_CPOFDM = hypergeom([1/2],[3/2 2],-(M_DopplerShift/F*pi).^2);
P_I_CPOFDM = 1-P_S_CPOFDM;
SIR_CPOFDM = 10*log10(P_S_CPOFDM./P_I_CPOFDM);


figure(21);
plot(M_Velocity_kmh,M_SIRoptimal(:,1),':','color', [0.5 0.5 0]);
hold on;
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,1)>5e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,1)>5e3,1),'color', [0.5 0.5 0]);
maxTemp = max(M_SIR007_Interpolated_True);
plot(M_Velocity_kmh(maxTemp>SIR_CPOFDM),maxTemp(maxTemp>SIR_CPOFDM),'black');  
plot(M_Velocity_kmh,SIR_CPOFDM,'red');
plot(M_Velocity_kmh,M_SIRoptimal(:,2),': blue');        
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,2)>5e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,2)>5e3,2),'blue');        
plot(M_Velocity_kmh,M_SIRoptimal(:,3),':','color',ColorFBMC);
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,3)>5e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,3)>5e3,3),'color',ColorFBMC);  

xlabel('Velocity [km/h]');
ylabel('Signal-to-Interference Ratio [dB]');
ylim([0 50]);
p1 = plot(nan,nan,': black');
p2 = plot(nan,nan,'black');
legend([p1 p2],{'Subcarrier Spacing, F=5kHz','Subcarrier Spacing, F>5kHz'},'Location','SouthWest');
    
if PlotAddionalResults
    % Optimal subcarrier spacing according to the rule of thumb   
    RMS_Doppler = 1/sqrt(2)*M_DopplerShift;
    RMS_Tau = 370e-9; % we have 44ns and not 46ns due to sampling 

    ApproxOptimalSubcarrierSpacing_PHYDYAS = sqrt(0.2745/0.328*RMS_Doppler/RMS_Tau);
    ApproxOptimalSubcarrierSpacing_Hermite = sqrt(0.2015/0.403*RMS_Doppler/RMS_Tau);

    figure(22);
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,1)>5e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,1)>5e3,1)/1e3,'red');
    hold on;
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,2)>5e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,2)>5e3,2)/1e3,'blue');        
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_Hermite/1e3,': blue');        
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,3)>5e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,3)>5e3,3)/1e3,'color',ColorFBMC);
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_PHYDYAS/1e3,':','color',ColorFBMC);        
    xlabel('Velocity [km/h], f_c=2.5GHz');
    ylabel('Optimal Subcarrier Spacing [kHz]');
    title('Vehicular A');
    p1 = plot(nan,nan,'black');
    p2 = plot(nan,nan,': black');
    legend([p1 p2],{'Exhaustive search','Approximation'},'Location','NorthWest');

    % SIR for approximate optimal subcarrier spacing
    rep1 = length(ApproxOptimalSubcarrierSpacing_PHYDYAS);
    rep2 = length(M_SubcarrierSpacingInterpolated);
    [~,index_PHYDYAS] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_PHYDYAS,[rep2,1])),[],1);
    [~,index_Hermite] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_Hermite,[rep2,1])),[],1);

    M_SIR_Interploated_PHYDYAS = M_SIR_Interpolated_ALL(:,:,3);
    M_SIR_Interploated_Hermite = M_SIR_Interpolated_ALL(:,:,2);
    for i_velocity = 1:length(M_Velocity_kmh)
        M_SIR_ApproxOpt_PHYDYAS(i_velocity) = M_SIR_Interploated_PHYDYAS(index_PHYDYAS(i_velocity),i_velocity);
        M_SIR_ApproxOpt_Hermite(i_velocity) = M_SIR_Interploated_Hermite(index_Hermite(i_velocity),i_velocity);

    end

    figure(23);plot(M_SIRoptimal(:,3)-M_SIR_ApproxOpt_PHYDYAS.');
    hold on;
    plot(M_SIRoptimal(:,2)-M_SIR_ApproxOpt_Hermite.');
    xlabel('Velocity [km/h], f_c=2.5GHz');
    ylabel('SIR Error in dB!!! between approx. and exhaustive search');
    title('Vehicular A');
    % SIR of theory vs simulations
    % plot(M_Velocity_kmh,[max(M_SIR007_Interpolated_True,[],1).' SIR_CPOFDM.']);
end


%% Plot TDL-A_30ns A
clear M_SIR_Interpolated_ALL M_SIRoptimal M_SIR_Interpolated M_SIR_15kHZ M_SubcarrierSpacingOptimal M_SIR007_Interpolated_True;
M_Velocity_kmh = [0 2 4 6 8 10 12 14 16 20:10:500];           
M_SubcarrierSpacingInterpolated = 1e3*[100:1:1000];
NrMethods = 4;
PowerDelayProfile = 'TDL-A_30ns';
Method{1} = 'OFDM';
PrototypeFilter{1} = '';        
Method{2} = 'FBMC-OQAM';
PrototypeFilter{2} = 'Hermite';
Method{3} = 'FBMC-OQAM';
PrototypeFilter{3} = 'PHYDYAS';
Method{4} = 'CP-OFDM-007LTE';
PrototypeFilter{4} = '';             

% Load Results and find optimum (using interpolation to obtain smoother curves)
for i_method = 1:NrMethods
    for i_velocity = 1:length(M_Velocity_kmh)   
        Velocity_kmh = M_Velocity_kmh(i_velocity);
        FileName = ['.\Results\' Method{i_method} '_' PrototypeFilter{i_method} '_' PowerDelayProfile int2str(Velocity_kmh) '.mat'];
        L1 = load(FileName);

        % Needed for interpolation: remove double subcarriers
        IndexSameSubcarrierSpacing = find(diff(L1.M_RealSubcarrierSpacing)==0);
        L1.M_RealSubcarrierSpacing(IndexSameSubcarrierSpacing) = [];
        L1.SIR_Theory(IndexSameSubcarrierSpacing) = [];

        % Interpolate to get smoother curves
        M_SIR_Interpolated(:,i_velocity) = interp1(L1.M_RealSubcarrierSpacing,L1.SIR_Theory,M_SubcarrierSpacingInterpolated,'spline',nan);

        if strcmp(Method{i_method},'CP-OFDM-007LTE') % To check if theoretical equations are correct (interpolation is wrong in case of CP-OFDM!)
            M_SIR007_Interpolated_True(:,i_velocity) = L1.SIR_Theory;
        end
    end
    M_SIR_Interpolated_ALL(:,:,i_method) = M_SIR_Interpolated;

    [SIRoptimal,IndexSIRoptimal] = max(M_SIR_Interpolated);
    SubcarrierSpacingOptimal = M_SubcarrierSpacingInterpolated(IndexSIRoptimal);

    M_SIR_100kHZ(:,i_method) = M_SIR_Interpolated(M_SubcarrierSpacingInterpolated==100e3,:);
    M_SIRoptimal(:,i_method) = SIRoptimal.';
    M_SubcarrierSpacingOptimal(:,i_method) = SubcarrierSpacingOptimal.';
end
M_SIR_100kHZ(:,4)=nan; % Values are wrong for CP-OFDM due to discontinuity
M_SIRoptimal(:,4)=nan;
M_SubcarrierSpacingOptimal(:,4)=nan;


% Calculate Theoretical SIR for OFDM (no ISI)
Kappa = 1/14;   % LTE efficiency loss
TCP =  289e-9;   % maximum delay
F = Kappa/TCP;  %optimal subcarrier spacing (no ISI)
M_DopplerShift = M_Velocity_kmh/3.6*60e9/2.998e8;
P_S_CPOFDM = hypergeom([1/2],[3/2 2],-(M_DopplerShift/F*pi).^2);
P_I_CPOFDM = 1-P_S_CPOFDM;
SIR_CPOFDM = 10*log10(P_S_CPOFDM./P_I_CPOFDM);


figure(31);
plot(M_Velocity_kmh,M_SIRoptimal(:,1),':','color', [0.5 0.5 0]);
hold on;
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,1)>100e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,1)>100e3,1),'color', [0.5 0.5 0]);
maxTemp = max(M_SIR007_Interpolated_True);
plot(M_Velocity_kmh(maxTemp>SIR_CPOFDM),maxTemp(maxTemp>SIR_CPOFDM),'black');  
plot(M_Velocity_kmh,SIR_CPOFDM,'red');
plot(M_Velocity_kmh,M_SIRoptimal(:,2),': blue');        
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,2)>100e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,2)>100e3,2),'blue');        
plot(M_Velocity_kmh,M_SIRoptimal(:,3),':','color',ColorFBMC);
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,3)>100e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,3)>100e3,3),'color',ColorFBMC);  


xlabel('Velocity [km/h]');
ylabel('Signal-to-Interference Ratio [dB]');
ylim([0 50]);
p1 = plot(nan,nan,': black');
p2 = plot(nan,nan,'black');
legend([p1 p2],{'Subcarrier Spacing, F=100kHz','Subcarrier Spacing, F>100kHz'},'Location','SouthWest');
 

if PlotAddionalResults
    % Optimal subcarrier spacing according to the rule of thumb   
    RMS_Doppler = 1/sqrt(2)*M_DopplerShift;
    RMS_Tau = 30e-9; % we have 44ns and not 46ns due to sampling 

    ApproxOptimalSubcarrierSpacing_PHYDYAS = sqrt(0.2745/0.328*RMS_Doppler/RMS_Tau);
    ApproxOptimalSubcarrierSpacing_Hermite = sqrt(0.2015/0.403*RMS_Doppler/RMS_Tau);

    figure(32);
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,1)>100e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,1)>100e3,1)/1e3,'red');
    hold on;
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,2)>100e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,2)>100e3,2)/1e3,'blue');        
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_Hermite/1e3,': blue');        
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,3)>100e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,3)>100e3,3)/1e3,'color',ColorFBMC);
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_PHYDYAS/1e3,':','color',ColorFBMC);        
    xlabel('Velocity [km/h], f_c=2.5GHz');
    ylabel('Optimal Subcarrier Spacing [kHz]');
    title('TDL-A');
    p1 = plot(nan,nan,'black');
    p2 = plot(nan,nan,': black');
    legend([p1 p2],{'Exhaustive search','Approximation'},'Location','NorthWest');

    % SIR for approximate optimal subcarrier spacing
    rep1 = length(ApproxOptimalSubcarrierSpacing_PHYDYAS);
    rep2 = length(M_SubcarrierSpacingInterpolated);
    [~,index_PHYDYAS] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_PHYDYAS,[rep2,1])),[],1);
    [~,index_Hermite] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_Hermite,[rep2,1])),[],1);

    M_SIR_Interploated_PHYDYAS = M_SIR_Interpolated_ALL(:,:,3);
    M_SIR_Interploated_Hermite = M_SIR_Interpolated_ALL(:,:,2);
    for i_velocity = 1:length(M_Velocity_kmh)
        M_SIR_ApproxOpt_PHYDYAS(i_velocity) = M_SIR_Interploated_PHYDYAS(index_PHYDYAS(i_velocity),i_velocity);
        M_SIR_ApproxOpt_Hermite(i_velocity) = M_SIR_Interploated_Hermite(index_Hermite(i_velocity),i_velocity);

    end

    figure(33);plot(M_SIRoptimal(:,3)-M_SIR_ApproxOpt_PHYDYAS.');
    hold on;
    plot(M_SIRoptimal(:,2)-M_SIR_ApproxOpt_Hermite.');
    xlabel('Velocity [km/h], f_c=60GHz$');
    ylabel('SIR Error in dB!!! between approx. and exhaustive search');
    title('TDL-A');
    % SIR of theory vs simulations
    % plot(M_Velocity_kmh,[max(M_SIR007_Interpolated_True,[],1).' SIR_CPOFDM.']);
end


%% Plot TDL-B_900ns A
clear M_SIR_ApproxOpt_Hermite M_SIR_Interpolated_ALL M_SIRoptimal M_SIR_Interpolated M_SIR_100kHZ M_SubcarrierSpacingOptimal M_SIR007_Interpolated_True;
M_Velocity_kmh = [0:10:500];        
M_SubcarrierSpacingInterpolated = 1e3*[100:1:1000];
NrMethods = 4;
PowerDelayProfile = 'TDL-B_900ns';
Method{1} = 'CP-OFDM';
PrototypeFilter{1} = '';        
Method{2} = 'FBMC-QAM';
PrototypeFilter{2} = 'Hermite';
Method{3} = 'FBMC-QAM';
PrototypeFilter{3} = 'PHYDYAS';
Method{4} = 'FBMC-OQAM';
PrototypeFilter{4} = 'Hermite';
     

% Load Results and find optimum (using interpolation to obtain smoother curves)
for i_method = 1:NrMethods
    for i_velocity = 1:length(M_Velocity_kmh)   
        Velocity_kmh = M_Velocity_kmh(i_velocity);
        FileName = ['.\Results\' Method{i_method} '_' PrototypeFilter{i_method} '_' PowerDelayProfile int2str(Velocity_kmh) '.mat'];
        L1 = load(FileName);

        % Needed for interpolation: remove double subcarriers
        IndexSameSubcarrierSpacing = find(diff(L1.M_RealSubcarrierSpacing)==0);
        L1.M_RealSubcarrierSpacing(IndexSameSubcarrierSpacing) = [];
        L1.SIR_Theory(IndexSameSubcarrierSpacing) = [];

        % Interpolate to get smoother curves
        M_SIR_Interpolated(:,i_velocity) = interp1(L1.M_RealSubcarrierSpacing,L1.SIR_Theory,M_SubcarrierSpacingInterpolated,'spline',nan);

        if strcmp(Method{i_method},'CP-OFDM') % To check if theoretical equations are correct (interpolation is wrong in case of CP-OFDM!)
            M_SIR2_Interpolated_True(:,i_velocity) = L1.SIR_Theory;
        end
    end
    M_SIR_Interpolated_ALL(:,:,i_method) = M_SIR_Interpolated;

    [SIRoptimal,IndexSIRoptimal] = max(M_SIR_Interpolated);
    SubcarrierSpacingOptimal = M_SubcarrierSpacingInterpolated(IndexSIRoptimal);

    M_SIR_100kHZ(:,i_method) = M_SIR_Interpolated(M_SubcarrierSpacingInterpolated==100e3,:);
    M_SIRoptimal(:,i_method) = SIRoptimal.';
    M_SubcarrierSpacingOptimal(:,i_method) = SubcarrierSpacingOptimal.';
end
M_SIR_100kHZ(:,1)=nan; % Values are wrong for CP-OFDM due to discontinuity
M_SIRoptimal(:,1)=nan;
M_SubcarrierSpacingOptimal(:,1)=nan;

% Calculate Theoretical SIR for OFDM (no ISI)
Kappa = 1;   % LTE efficiency loss
TCP =  4305e-9;   % maximum delay
F = Kappa/TCP;  %optimal subcarrier spacing (no ISI)
M_DopplerShift = M_Velocity_kmh/3.6*60e9/2.998e8;
P_S_CPOFDM = hypergeom([1/2],[3/2 2],-(M_DopplerShift/F*pi).^2);
P_I_CPOFDM = 1-P_S_CPOFDM;
SIR_CPOFDM = 10*log10(P_S_CPOFDM./P_I_CPOFDM);


figure(41);
maxTemp = max(M_SIR2_Interpolated_True);
plot(M_Velocity_kmh(maxTemp>SIR_CPOFDM),maxTemp(maxTemp>SIR_CPOFDM),'black');  
hold on;
plot(M_Velocity_kmh,SIR_CPOFDM,'red');
plot(M_Velocity_kmh,M_SIRoptimal(:,2),': blue');        
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,2)>100e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,2)>100e3,2),'blue');        
plot(M_Velocity_kmh,M_SIRoptimal(:,3),':','color',ColorFBMC);
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,3)>100e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,3)>100e3,3),'color',ColorFBMC);  
plot(M_Velocity_kmh,M_SIRoptimal(:,4),': blue');        
plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,4)>100e3),M_SIRoptimal(M_SubcarrierSpacingOptimal(:,4)>100e3,4),'blue');        


xlabel('Velocity [km/h]');
ylabel('Signal-to-Interference Ratio [dB]');
ylim([0 50]);
p1 = plot(nan,nan,': black');
p2 = plot(nan,nan,'black');
legend([p1 p2],{'Subcarrier Spacing, F=100kHz','Subcarrier Spacing, F>100kHz'},'Location','SouthWest');


if PlotAddionalResults
    % Optimal subcarrier spacing according to the rule of thumb   
    RMS_Doppler = 1/sqrt(2)*M_DopplerShift;
    RMS_Tau = 900e-9; % we have 44ns and not 46ns due to sampling 

    ApproxOptimalSubcarrierSpacing_Hermite = sqrt(0.2015/0.403*RMS_Doppler/RMS_Tau);
    ApproxOptimalSubcarrierSpacing_PHYDYAS_TF2 = sqrt(0.2745/0.328*RMS_Doppler/RMS_Tau*2); % in FBMC-QAM we have F=2/T_0 => compared to OQAM (F=1/T_0) we need an additional factor of two
    ApproxOptimalSubcarrierSpacing_Hermite_TF2 = sqrt(0.2015/0.403*RMS_Doppler/RMS_Tau*2);  % in FBMC-QAM we have F=2/T_0 => compared to OQAM (F=1/T_0) we need an additional factor of two

    figure(42);
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,2)>100e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,2)>100e3,2)/1e3,'blue');        
    hold on;
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_Hermite_TF2/1e3,': blue');        
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,3)>100e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,3)>100e3,3)/1e3,'color',ColorFBMC);
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_PHYDYAS_TF2/1e3,':','color',ColorFBMC); 
    plot(M_Velocity_kmh(M_SubcarrierSpacingOptimal(:,4)>100e3),M_SubcarrierSpacingOptimal(M_SubcarrierSpacingOptimal(:,4)>100e3,4)/1e3,'blue');
    plot(M_Velocity_kmh,ApproxOptimalSubcarrierSpacing_Hermite/1e3,':','color',ColorFBMC);     
    
    xlabel('Velocity [km/h], f_c=2.5GHz');
    ylabel('Optimal Subcarrier Spacing [kHz]');
    title('TDL-B 900ns');
    p1 = plot(nan,nan,'black');
    p2 = plot(nan,nan,': black');
    legend([p1 p2],{'Exhaustive search','Approximation'},'Location','NorthWest');

    % SIR for approximate optimal subcarrier spacing
    rep1 = length(ApproxOptimalSubcarrierSpacing_PHYDYAS_TF2);
    rep2 = length(M_SubcarrierSpacingInterpolated);
    [~,index_PHYDYAS_TF2] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_PHYDYAS_TF2,[rep2,1])),[],1);
    [~,index_Hermite_TF2] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_Hermite_TF2,[rep2,1])),[],1);
    [~,index_Hermite] = min(abs(repmat(M_SubcarrierSpacingInterpolated.',[1 rep1])-repmat(ApproxOptimalSubcarrierSpacing_Hermite,[rep2,1])),[],1);

 
    M_SIR_Interploated_PHYDYAS_TF2 = M_SIR_Interpolated_ALL(:,:,3);
    M_SIR_Interploated_Hermite_TF2 = M_SIR_Interpolated_ALL(:,:,2);
    M_SIR_Interploated_Hermite = M_SIR_Interpolated_ALL(:,:,4);
    for i_velocity = 1:length(M_Velocity_kmh)
        M_SIR_ApproxOpt_PHYDYAS_TF2(i_velocity) = M_SIR_Interploated_PHYDYAS_TF2(index_PHYDYAS_TF2(i_velocity),i_velocity);
        M_SIR_ApproxOpt_Hermite_TF2(i_velocity) = M_SIR_Interploated_Hermite_TF2(index_Hermite_TF2(i_velocity),i_velocity);
        M_SIR_ApproxOpt_Hermite(i_velocity) = M_SIR_Interploated_Hermite(index_Hermite(i_velocity),i_velocity);

    end

    figure(43);plot(M_SIRoptimal(:,3)-M_SIR_ApproxOpt_PHYDYAS_TF2.');
    hold on;
    plot(M_SIRoptimal(:,2)-M_SIR_ApproxOpt_Hermite_TF2.');
    plot(M_SIRoptimal(:,4)-M_SIR_ApproxOpt_Hermite.');
    xlabel('Velocity [km/h], f_c=60GHz');
    ylabel('SIR Error in dB!!! between approx. and exhaustive search');
    title('TDL-B 900ns');
    % SIR of theory vs simulations
    % plot(M_Velocity_kmh,[max(M_SIR2_Interpolated_True,[],1).' SIR_CPOFDM.']);
end



