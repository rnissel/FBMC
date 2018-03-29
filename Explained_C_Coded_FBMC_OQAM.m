% =========================================================================   
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.nt.tuwien.ac.at
% =========================================================================      
% This script shows how to find the precoding matrix which allows QAM
% transmission in FBMC-OQAM at full rate. Furthermore, it illustrates the
% concept of time/frequency spreading.

clear; close all;

%% Parameters
FrequencySpreading = true;                                      % true = frequency spreading; false = time spreading
L = 2^4;                                                        % Number of subcarriers (must be power of 2 if FrequencySpreading==true)
K = 2^5;                                                        % Number of FBMC symbols in time (must be power of 2 if FrequencySpreading==false)

%% FBMC object
FBMC = Modulation.FBMC(...
    L,...                                                       % Number subcarriers
    K,...                                                       % Number FBMC symbols (determines the frequency resolution!)
    15e3,...                                                    % Subcarrier spacing (Hz)
    15e3*2*L,...                                                % Sampling rate (Samples/s)
    0,...                                                       % Intermediate frequency first subcarrier (Hz)
    false,...                                                   % Transmit real valued signal
    'Hermite-OQAM',...                                          % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    8, ...                                                      % Overlapping factor (also determines oversampling in the frequency domain)
    0, ...                                                      % Initial phase shift
    true ...                                                    % Polyphase implementation
    );

% Transmit Matrix G, See also "A_FBMCOQAM_Explained.m"
G = FBMC.GetTXMatrix;
G = G/sqrt(G(:,1)'*G(:,1)); % Normalize

% Find coding matrix C which allows QAM transmissions in FBMC-OQAM
if FrequencySpreading 
    % see below (28)
    H   = fwht(eye(L),L,'sequency')*sqrt(L); % Fast Walsh-Hadamard Transform
    C0  = H(:,1:2:end);
    C   = kron(eye(K),C0);
else % Time spreading
    C   = FBMC.GetPrecodingMatrixForQAMinOQAM(1,1);
end
  
% Check orthogonality
figure(1);
imagesc(10*log10(abs(C'*G'*G*C).^2),[-50 0]);
title('Check orthogonality: C''G''GC (should be diagonal)');


% Illustration of the spreading approach. We plot, for each step, one 
% column of the reshaped coding matrix C, showing over which time-frequency
% positions the data symbol is spread
minminC = min(min(C));
maxmaxC = max(max(C));
if FrequencySpreading
    for i_lk = 1:L/2*K
            figure(2);
            imagesc(reshape(C(:,i_lk),[L,K]),[minminC maxmaxC]);
            xlabel('Time-position, k');
            ylabel('Frequency-position, l');
            title('Spreading vector');
            pause(0.1);
    end
else
    IndexTemp = reshape(1:L*K/2,L,K/2);
    for i_l=1:L
        for i_k = 1:K/2 
            figure(2);
            imagesc(reshape(C(:,IndexTemp(i_l,i_k)),[L,K]),[minminC maxmaxC]);
            xlabel('Time-position, k');
            ylabel('Frequency-position, l');
            title('Spreading vector');            
            pause(0.1);            
        end
    end       
end

