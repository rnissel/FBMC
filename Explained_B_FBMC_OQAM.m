% =========================================================================   
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =========================================================================      
% This script describes a back-to-back FBMC-OQAM transmission, that is,
% without channel. We consider two possible implementations
%       1) Based on a transmit matrix G
%       2) Efficient IFFT implementation
% The output compares the received symbols to the transmitted symbols.
% Furthermore, we plot the transmit power over time and the power spectral
% density. 

clear; close all;

%% Parameters
fs  = 1024*15e3;                 % Sampling frequency
F   = 15e3;                      % Subcarrier spacing (frequency spacing)
O   = 4;                         % Overlapping factor
L   = 4;                         % Number of subcarriers
K   = 3;                         % Number of FBMC symbols in time (each consists of L subcarriers)

%% Dependent Parameters
dt  = 1/fs;                      % Time between samples


% We choose the Hermite prototype filter because it is more flexible in terms of overlapping factor. For other pulses, see "A_PrototypeFilters.m"
% See (10) and (11) for the definition of the Hermite prototype filter
p_Hermite = @(t,T0) ((t<=(O*T0/2))&(t>-(O*T0/2))).* ...
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
                     
% For OQAM
T0  = 1/F;          % Time-scaling parameter
T   = T0/2;         % Time spacing


% Number of total samples, round due to numerical inaccuracies
N = round((O*T0+T*(K-1))*fs); 
t = (0:N-1)*dt-O*T0/2;

%% Transmit Matrix G, see (2), (18) and (22)
for l=0:L-1
    for k=0:K-1
        G(:,l+1,k+1)=p_Hermite(t-k*T,T0).*exp(1j*2*pi*l*F*(t-k*T)).*exp(1j*pi/2*(l+k))*sqrt(dt);
    end
end
G=G(:,:);


%% IFFT Implementation: Transmitter
N_FFT   = round(fs/F);                                          % FFT size, round to avoid numerical inaccuracies
W       = dftmtx(N_FFT);                                        % DFT matrix
p       = sqrt(dt)*(p_Hermite((0:N_FFT*O-1)*dt-O*T0/2,T0)).';   % Sampled prototype filter.

% Generate random transmit symbols
x_2D    = (randi(2,L,K,1)-1.5)*2; % Generate random data symbols of size L \times K
x       = x_2D(:);                % Vectorized data symbols

% Calculate the samples for each time position k separately and map them to the corresponding time-position 
s_IFFT           = zeros(N,K);
s_MatrixNotation = zeros(N,K);
for k = 0:K-1
    TimePositionIndex = (1:O*N_FFT)+k*N_FFT/2;
    % Matrix representation as in the paper, see (32)
    s_MatrixNotation(TimePositionIndex,k+1) = p.*(...
                          kron(...
                                ones(O,1),...
                                W'*[x_2D(:,k+1).*1j.^((0:L-1)'+k);zeros(N_FFT-L,1)]...
                              )...
                        );
    % Equivalent calculation. Note that we have to scale the IFFT by N_FFT to be consistent with the matrix notation
    s_IFFT(TimePositionIndex,k+1) = p.*repmat(ifft(x_2D(:,k+1).*1j.^((0:L-1)'+k),N_FFT),O,1)*(N_FFT);
end
% Add the overlapping symbols in time
s_MatrixNotation    = sum(s_MatrixNotation,2);
s_IFFT              = sum(s_IFFT,2); % should be the same as "MatrixNotation"


%% FFT Implementation: Receiver (just the reversed of the transmitter)
for k = 0:K-1
    TimePositionIndex = (1:O*N_FFT)+k*N_FFT/2;  
    y_temp = fft(sum(reshape(p.*s_IFFT(TimePositionIndex),N_FFT,O),2));
    y_IFFT(:,k+1) = y_temp(1:L).*1j.^(-((0:L-1)'+k));
end


%% Check results: real part is the same but we observe imaginary interference
disp('  Received symbols  | Received symbols | Transmitted symbols');
disp('     y (using G)    |  y (using IFFT)  |        x     ');
disp([      G'*G*x             y_IFFT(:)              x       ]);


%% Plot Transmit Power
Ps = sum(G.*conj(G),2); % same as: diag(G*G');
figure();
plot(t,Ps);
xlabel('Time (s)');
ylabel('Transmit Power (not normalized)');


%% Plot Power Spectral Density
PSD = sum(abs(fft(G)).^2,2);
PSD = PSD/max(PSD);
f = (0:N-1)*fs/N;
figure();
plot(f,10*log10(abs(PSD)))
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');



