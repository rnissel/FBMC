% =========================================================================   
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =========================================================================      
% This script plots the time and frequency response of different prototype
% filters:
%       1) Rectangular
%       2) PHYDYAS, Overlapping factor 4 (but 8 is also implemented)
%       3) Hermite
%       4) RRC in time

clear; close all;

fs                      = 15.36e6;          % Sampling Rate: 1024*15kHz
F                       = 15e3;             % Subcarrier Spacing: 15kHz

dt = 1/fs;

%% Define Prototype Filters! Orthogonal for T=T0 and F=2/T0!
% PHYDYAS pulse for overlapping factor 4
p_PHYDYAS_O4 = @(t,T0) ((t<=(4*T0/2))&(t>-4*T0/2)).*(1+2*(...
    0.97195983 * cos(2*pi*1/4*t/T0) + ...
    sqrt(2)/2  * cos(2*pi*2/4*t/T0) + ...
    0.23514695 * cos(2*pi*3/4*t/T0) ...
    ))/sqrt(4^2*T0);

% PHYDYAS pulse for overlapping factor 8
p_PHYDYAS_O8 = @(t,T0) ((t<=(8*T0/2))&(t>-8*T0/2)).*(1+2*(...
    0.99932588 * cos(2*pi*1/8*t/T0) + ...
    0.98203168 * cos(2*pi*2/8*t/T0) + ...
    0.89425129 * cos(2*pi*3/8*t/T0) + ...
    sqrt(2)/2  * cos(2*pi*4/8*t/T0) + ...
    0.44756522 * cos(2*pi*5/8*t/T0) + ...
    0.18871614 * cos(2*pi*6/8*t/T0) + ...
    0.03671221 * cos(2*pi*7/8*t/T0) ...
    ))/sqrt(8^2*T0);

% Root raised cosine filter in time (low latency)
p_timeRRC = @(t,T0) ((t<=(T0/2))&(t>-T0/2)).*sqrt(1+(...
      cos(pi*1*2*t/T0)  ...
     ))/sqrt(T0);
 
% Hermite prototype filter for overlapping factor 8
O  = 8; 
p_Hermite = @(t,T0) ((t<=(O*T0/2))&(t>-O*T0/2)).* ...
    1/sqrt(T0).*exp(-pi*(t./(T0/sqrt(2))).^2) .* (...
    1.412692577 + ...
    -3.0145e-3 .*...
                ((12+(-48).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+16.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4 ) )+ ...
    -8.8041e-6 .*...
                (1680+(-13440).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+13440.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4+(-3584).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^6+256.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^8 )+ ...
    -2.2611e-9  .*... 
                 (665280+(-7983360).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+13305600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4+(-7096320).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^6+1520640.* ...
                    (sqrt(2*pi)*(t./(T0/sqrt(2)))).^8+(-135168).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^10+4096.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^12 )+ ...
    -4.4570e-15 .*... 
                 (518918400+(-8302694400).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+19372953600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4+(-15498362880).* ...
                   (sqrt(2*pi)*(t./(T0/sqrt(2)))).^6+5535129600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^8+(-984023040).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^10+89456640.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^12+( ...
                   -3932160).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^14+65536.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^16 )+ ...
     1.8633e-16 .*...
                 (670442572800+(-13408851456000).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^2+40226554368000.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^4+( ...
                   -42908324659200).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^6+21454162329600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^8+(-5721109954560).* ...
                   (sqrt(2*pi)*(t./(T0/sqrt(2)))).^10+866834841600.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^12+(-76205260800).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^14+3810263040.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^16+ ...
                   (-99614720).*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^18+1048576.*(sqrt(2*pi)*(t./(T0/sqrt(2)))).^20 ));
                
% Reference rectangular prototype filter, orthogonal for T=T0, F=1/T0!
p_Rectangular = @(t,T0) (1/sqrt(T0).*((t<=(T0/2))&(t>-T0/2)));


%% Plot Filters
T0  = 1/F;
t   = (-200*T0/2:dt:(200*T0/2)).'; t(end)=[]; % so that fft becomes real

% Time domain
p_PHYDYAS_O4_Samples    = p_PHYDYAS_O4(t,T0);
p_PHYDYAS_O8_Samples    = p_PHYDYAS_O8(t,T0);
p_timeRRC_Samples       = p_timeRRC(t,T0);
p_Hermite_Samples       = p_Hermite(t,T0);
p_Rectangular_Samples   = p_Rectangular(t,T0);

% Plot the prototype filters in time
figure();
plot(t/T0,p_PHYDYAS_O8_Samples*sqrt(T0),'red');
hold on;
plot(t/T0,p_timeRRC_Samples*sqrt(T0),'magenta');
plot(t/T0,p_Hermite_Samples*sqrt(T0),'blue');
plot(t/T0,p_Rectangular_Samples*sqrt(T0),'black');
xlim([-2 2]);
ylabel('p(t)');
xlabel('Normalized Time, t/T0');
legend({'PHYDYAS','RRC','Hermite','Rectangular'});

% Frequency domain
df  = 1/(length(t)*dt);
f   = (-length(t)/2:(length(t)/2-1))*df; 

P_PHYDYAS_O8_FFT    = circshift(abs(fft(p_PHYDYAS_O8_Samples)),[-length(t)/2,0]);
P_PHYDYAS_O4_FFT    = circshift(abs(fft(p_PHYDYAS_O4_Samples)),[-length(t)/2,0]);
P_timeRRC_FFT       = circshift(abs(fft(p_timeRRC_Samples)),[-length(t)/2,0]);
P_Hermite_FFT       = circshift(abs(fft(p_Hermite_Samples)),[-length(t)/2,0]);
P_Rectangular_FFT   = circshift(abs(fft(p_Rectangular_Samples)),[-length(t)/2,0]);

Normalize           = max(max([P_PHYDYAS_O4_FFT P_timeRRC_FFT P_Hermite_FFT P_Rectangular_FFT]));

% Plot the prototype filters in frequency
figure();
plot(f*T0,20*log10(P_PHYDYAS_O4_FFT/Normalize),'red');
hold on;
plot(f*T0,20*log10(P_timeRRC_FFT/Normalize),'magenta');
plot(f*T0,20*log10(P_Hermite_FFT/Normalize),'blue');
plot(f*T0,20*log10(P_Rectangular_FFT/Normalize),'black');
ylim([-100 0]);
xlim([-10 10]);
ylabel('P(f)');
xlabel('Normalized Frequency, f T0');
legend({'PHYDYAS','RRC','Hermite','Rectangular'});

