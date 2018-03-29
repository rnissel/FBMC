classdef FOFDM < handle
    % =====================================================================        
    % This MATLAB class represents an implementation of FOFDM (filtered OFDM). 
    % The modulation parameters are initialized by the class contructor. 
    % The modulation of data symbols x and the demodulation of the
    % received samples r is then performed by the methods ".Modulation(x)"
    % and ".Demodulation(r)".
    % Additionally, there exist some other useful methods, such as,
    % ".PlotPowerSpectralDensity" or ".GetTXMatrix".
    % =====================================================================    
    % Ronald Nissel, rnissel@nt.tuwien.ac.at
    % (c) 2017 by Institute of Telecommunications, TU Wien
    % www.nt.tuwien.ac.at
    % =====================================================================    
    
    properties (SetAccess = private)
        Nr                  % for dimensionless parameters
        PHY                 % for parameters with physical interpretation
        Implementation      % implmentation relevent parameters
    end
    
    methods
        function obj = FOFDM(varargin)
            % Initialize parameters, set default values
            if numel(varargin) == 11
                obj.Nr.Subcarriers               = varargin{1};     % Number of subcarriers
                obj.Nr.MCSymbols                 = varargin{2};     % Number FOFDM symbols in time
                obj.PHY.SubcarrierSpacing        = varargin{3};     % Subcarrier spacing (Hz)
                obj.PHY.SamplingRate             = varargin{4};     % Sampling rate (Samples/s)
                obj.PHY.IntermediateFrequency    = varargin{5};     % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
                obj.PHY.TransmitRealSignal       = varargin{6};     % Transmit real valued signal (sampling theorem must be fulfilled!)
                obj.PHY.CyclicPrefixLength       = varargin{7};     % Length of the cyclic prefix (s)
                obj.PHY.ZeroGuardTimeLength      = varargin{8};     % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
                % FOFDM
                obj.PHY.FilterLengthTX           = varargin{9};     % Length of the transmit filter (s)
                obj.PHY.FilterLengthRX           = varargin{10};    % Length of the receive filter (s)  
                obj.PHY.CyclicPrefixLengthFilter = varargin{11};    % Length of the additional cyclic prefix (s).  Needed to combat ISI and ICI due to the filtering. However, some small ICI and ISI is perfectly fine.
            elseif numel(varargin) == 0
                % LTE like values (can later be changed using OFDM.Set...)
                obj.Nr.Subcarriers               = 24;
                obj.Nr.MCSymbols                 = 14;
                obj.PHY.SubcarrierSpacing        = 15e3;
                obj.PHY.SamplingRate             = 15e3*24*14;
                obj.PHY.IntermediateFrequency    = 0;
                obj.PHY.TransmitRealSignal       = false;
                obj.PHY.CyclicPrefixLength       = 0;
                obj.PHY.ZeroGuardTimeLength      = 0;
                % FOFDM
                obj.PHY.FilterLengthTX           = 1.5/(14*15e3);        % Chose so that no ISI/ICI. However, some small ICI/ISI is usually perfectly fine.  
                obj.PHY.FilterLengthRX           = 1.5/(14*15e3);        % Chose so that no ISI/ICI. However, some small ICI/ISI is usually perfectly fine. 
                obj.PHY.CyclicPrefixLengthFilter = 1/(14*15e3);                  
            else
                error('Number of input variables must be either 0 (default values) or 11');
            end
            
            % calculate and set all dependent parameters
            obj.SetDependentParameters();            
            
        end
        
        function SetDependentParameters(obj)
            % method that sets all parameters which are dependent on other
            % parameters
            
            if mod(round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)*10^5)/10^5,1)~=0
                obj.PHY.SubcarrierSpacing=obj.PHY.SamplingRate/(round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)));
                disp('Sampling rate must be a multiple of the subcarrier spacing!');
                disp(['Therefore, the subcarrier spacing is set to: ' int2str(obj.PHY.SubcarrierSpacing) 'Hz']);
            end
            
            if mod(round(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing*10^5)/10^5,1)~=0
                obj.PHY.IntermediateFrequency = round(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing)*obj.PHY.SubcarrierSpacing;
                disp('The intermediate frequency must be a multiple of the subcarrier spacing!');
                disp(['Therefore, the intermediate frequency is set to ' int2str(obj.PHY.IntermediateFrequency) 'Hz']);
            end
            
            if (obj.PHY.SamplingRate<obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing)
                error('Sampling theorem is not fullfilled: sampling rate must be higher than the number of subcarriers times subcarrier spacing');
            end
            
            if abs(mod(round(obj.PHY.CyclicPrefixLength*obj.PHY.SamplingRate*10^5)/10^5,1))~=0
                obj.PHY.CyclicPrefixLength=round(obj.PHY.CyclicPrefixLength*obj.PHY.SamplingRate)/obj.PHY.SamplingRate;
                disp('The length of the cyclic prefix times the sampling rate must be an integer!');
                disp(['Therefore, the cyclic prefix length is set to ' num2str(obj.PHY.CyclicPrefixLength) 's']);
            end
            
            obj.Implementation.CyclicPrefix             = round( obj.PHY.CyclicPrefixLength * obj.PHY.SamplingRate );                                       % number of samples for the CP
            obj.Implementation.ZeroGuardSamples         = round( obj.PHY.ZeroGuardTimeLength * obj.PHY.SamplingRate );                                      % number of samples for the zero guard
%             obj.Implementation.TimeSpacing              = round( obj.PHY.SamplingRate / obj.PHY.SubcarrierSpacing ) + obj.Implementation.CyclicPrefix;      % number of samples per sybol including CP
            obj.Implementation.FFTSize                  = round( obj.PHY.SamplingRate / obj.PHY.SubcarrierSpacing );
            obj.Implementation.IntermediateFrequency    = round( obj.PHY.IntermediateFrequency / obj.PHY.SubcarrierSpacing );
            obj.Implementation.NormalizationFactor      = sqrt( obj.PHY.SamplingRate^2 / obj.PHY.SubcarrierSpacing^2 / obj.Nr.Subcarriers);                 % Normalization factor so that power = 1
            obj.PHY.dt                                  = 1 / obj.PHY.SamplingRate;                                                                         % inter sample time
%             obj.PHY.TimeSpacing                         = obj.Implementation.TimeSpacing * obj.PHY.dt;                                                      % symbol duration including CP
%             obj.Nr.SamplesTotal                         = (obj.Nr.MCSymbols*obj.Implementation.TimeSpacing)+2*obj.Implementation.ZeroGuardSamples;          % total number of samples for all symbols including zero guards
        

            % FOFDM
            obj.Implementation.FilterLengthTX           = round(obj.PHY.FilterLengthTX *obj.PHY.SamplingRate)+3; % 3 because first and last filter elements are zero due to hann windowing. Futhermore we need one sample for the dirac pulse ;)
            obj.Implementation.FilterLengthRX           = round(obj.PHY.FilterLengthRX *obj.PHY.SamplingRate)+3;
            obj.Implementation.CyclicPrefixLengthFilter = round(obj.PHY.CyclicPrefixLengthFilter*obj.PHY.SamplingRate);
            
            obj.PHY.FilterLengthTX           = (obj.Implementation.FilterLengthTX-3) *obj.PHY.dt ;
            obj.PHY.FilterLengthRX           = (obj.Implementation.FilterLengthRX-3) *obj.PHY.dt;
            obj.PHY.CyclicPrefixLengthFilter = obj.Implementation.CyclicPrefixLengthFilter *obj.PHY.dt;
                        


            % FOFDM: sinc-filter + Hann time-window. 
            GuardSubcarriers = 0;
            FilterFrequency = circshift([ones(obj.Nr.Subcarriers+GuardSubcarriers,1);...
                 zeros(obj.Implementation.FFTSize-obj.Nr.Subcarriers-2*GuardSubcarriers,1);ones(GuardSubcarriers,1)],...
                 [obj.Implementation.IntermediateFrequency,0]);
            FilterLong = ifft(FilterFrequency);   
            obj.Implementation.FilterImpulseResponseTX = hann(obj.Implementation.FilterLengthTX).*[FilterLong(end-floor(obj.Implementation.FilterLengthTX/2)+1:end);FilterLong(1:ceil(obj.Implementation.FilterLengthTX/2))];
            obj.Implementation.FilterImpulseResponseRX = hann(obj.Implementation.FilterLengthRX).*[FilterLong(end-floor(obj.Implementation.FilterLengthRX/2)+1:end);FilterLong(1:ceil(obj.Implementation.FilterLengthRX/2))];


            obj.Implementation.TimeSpacing              = round( obj.PHY.SamplingRate / obj.PHY.SubcarrierSpacing ) + obj.Implementation.CyclicPrefix + obj.Implementation.CyclicPrefixLengthFilter;      % number of samples per sybol including CP
            obj.PHY.TimeSpacing                         = obj.Implementation.TimeSpacing * obj.PHY.dt;                                                      % symbol duration including CP
  
            obj.Nr.SamplesTotal = (obj.Nr.MCSymbols*obj.Implementation.TimeSpacing)+2*obj.Implementation.ZeroGuardSamples+obj.Implementation.FilterLengthTX-1;             
          
            FilterTransferFunctionTX = fft([obj.Implementation.FilterImpulseResponseTX;zeros(obj.Implementation.FFTSize-obj.Implementation.FilterLengthTX,1)]);
            FilterTransferFunctionRX = fft([obj.Implementation.FilterImpulseResponseRX;zeros(obj.Implementation.FFTSize-obj.Implementation.FilterLengthRX,1)]);
  
            obj.Implementation.FilterPreEqualizer = 1./FilterTransferFunctionTX(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers));
            obj.Implementation.FilterEqualizerRX  = 1./FilterTransferFunctionRX(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers));
            
            % Change normalization factor in order to account for Filter!
            TXMatrix1OFDMSymbolTemp = ifft(diag(circshift([obj.Implementation.FilterPreEqualizer;zeros(obj.Implementation.FFTSize-obj.Nr.Subcarriers,1)],[obj.Implementation.IntermediateFrequency,0])))*obj.Implementation.NormalizationFactor;
            TXMatrix1OFDMSymbolTemp = conv2([TXMatrix1OFDMSymbolTemp(end-obj.Implementation.CyclicPrefix-obj.Implementation.CyclicPrefixLengthFilter+1:end,:);TXMatrix1OFDMSymbolTemp],obj.Implementation.FilterImpulseResponseTX);

            obj.Implementation.FilterPowerNormalization = 1/(sum(sum(abs(TXMatrix1OFDMSymbolTemp).^2,2))/(obj.Implementation.FFTSize+obj.Implementation.CyclicPrefix+obj.Implementation.CyclicPrefixLengthFilter));
            obj.Implementation.NormalizationFactor = obj.Implementation.NormalizationFactor*sqrt(obj.Implementation.FilterPowerNormalization);    
        end
        
        % Set Functions
        function SetNrSubcarriers(varargin) 
            % set the number of subcarriers
            
            obj = varargin{1};
            % set specific property
            obj.Nr.Subcarriers = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetNrMCSymbols(varargin)
            % set the number of symbols
            
            obj = varargin{1};
            % set specific property
            obj.Nr.MCSymbols = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetSubcarrierSpacing(varargin)
            % set the subcarrier spacing
            
            obj = varargin{1};
            % set specific property
            obj.PHY.SubcarrierSpacing = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetSamplingRate(varargin)
            % set the sampling rate
            
            obj = varargin{1};
            % set specific property
            obj.PHY.SamplingRate = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetIntermediateFrequency(varargin)
            % set the intermediate frequency
            
            obj = varargin{1};
            % set specific property
            obj.PHY.IntermediateFrequency = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        function SetTransmitRealSignal(varargin)
            % set real transmit signal indicator
            
            obj = varargin{1};
            % set specific property
            obj.PHY.TransmitRealSignal = varargin{2};
            % recalculate dependent parameters
            obj.SetDependentParameters;
        end
        
        
        % Modulation and Demodulation
        function TransmitSignal = Modulation(obj,DataSymbols)
            % modulate the data symbols. The input argument is a matrix
            % of size "Number of subcarriers" \times "Number of FOFDM symbols (time)"
            % which represents the transmit data symbols
            
            DataSymbolsTemp = zeros(obj.Implementation.FFTSize,obj.Nr.MCSymbols);
            DataSymbolsTemp(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers),:) =  bsxfun(@times,DataSymbols,obj.Implementation.FilterPreEqualizer)*obj.Implementation.NormalizationFactor;
            if obj.PHY.TransmitRealSignal
                DataSymbolsTemp = (DataSymbolsTemp+conj(DataSymbolsTemp([1 end:-1:2],:)))/sqrt(2);
            end
            TransmitSignalNoCP = ifft(DataSymbolsTemp);
            TransmitSignal = [zeros(obj.Implementation.ZeroGuardSamples,1);reshape([TransmitSignalNoCP(end-obj.Implementation.CyclicPrefix-obj.Implementation.CyclicPrefixLengthFilter+1:end,:);TransmitSignalNoCP],obj.Implementation.TimeSpacing*obj.Nr.MCSymbols,1);zeros(obj.Implementation.ZeroGuardSamples,1)];
 
            TransmitSignal = conv(TransmitSignal,obj.Implementation.FilterImpulseResponseTX);

        end
        
        function ReceivedSymbols = Demodulation(obj, ReceivedSignal)
            % demodulates the received time signal and returns a matrix of 
            % size "Number of subcarriers" \times "Number of FOFDM symbols"
            % which represents the received symbols after demodulation but
            % before channel equalization
            
            ReceivedSignalAfterFilter = conv(ReceivedSignal,obj.Implementation.FilterImpulseResponseRX);
            
            % We only have a CP and not a CP+CS => if the CP is
            % not sufficiently long, we should timeshift the signal so that
            % interference occurse at both sides and not just at one. This
            % must later be compensated by circshift in order to get the
            % correct phase.
            FilterLengthMinusFilterCPHalf = max([floor((obj.Implementation.FilterLengthTX+obj.Implementation.FilterLengthRX-obj.Implementation.CyclicPrefixLengthFilter)/2)-1 2]); % 2 because two times zero elements for tx and rx filter
            
            ReceivedSignal_reshape = reshape(ReceivedSignalAfterFilter((obj.Implementation.ZeroGuardSamples+1+FilterLengthMinusFilterCPHalf):(end-obj.Implementation.ZeroGuardSamples-obj.Implementation.FilterLengthTX-obj.Implementation.FilterLengthRX+2+FilterLengthMinusFilterCPHalf)),obj.Implementation.TimeSpacing,obj.Nr.MCSymbols);
            ReceivedSignal_reshape_BeforeFFT = ReceivedSignal_reshape(obj.Implementation.CyclicPrefix+obj.Implementation.CyclicPrefixLengthFilter+1:end,:);
            ReceivedSymbolsTemp = fft(circshift(ReceivedSignal_reshape_BeforeFFT,[FilterLengthMinusFilterCPHalf 0]));
            if obj.PHY.TransmitRealSignal
                ReceivedSymbolsTemp = ReceivedSymbolsTemp*sqrt(2);
            end
            % Compensate also the receive filter
            ReceivedSymbols = ReceivedSymbolsTemp(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers),:)/(obj.Implementation.NormalizationFactor).*repmat(obj.Implementation.FilterEqualizerRX,1,obj.Nr.MCSymbols);  
        end
        
        % Matrix Description
        function TXMatrix = GetTXMatrix(obj)
            % returns a matrix G so that s=G*x(:) is the same as 
            % s=obj.Modulation(x)
            
            if obj.PHY.TransmitRealSignal
                % Does not work because conjungate complex is not a linear operation for complex symols. For FBMC it workes because real symbols are transmitted.
                error('GetTXMatrix is not supported for PHY.TransmitRealSignal == true!')
            end
            TXMatrix=zeros(obj.Nr.SamplesTotal,obj.Nr.Subcarriers*obj.Nr.MCSymbols);
            TXMatrixTemp=zeros(obj.Nr.SamplesTotal,obj.Nr.Subcarriers);
            x = zeros(obj.Nr.Subcarriers, obj.Nr.MCSymbols);
            for i_l= 1:obj.Nr.Subcarriers;
                x(i_l)=1;
                TXMatrixTemp(:,i_l) = obj.Modulation(x);
                x(i_l)=0;
            end
            for i_k=1:obj.Nr.MCSymbols
                TXMatrix(:,(1:obj.Nr.Subcarriers)+(i_k-1)*obj.Nr.Subcarriers)=circshift(TXMatrixTemp,[(i_k-1)*obj.Implementation.TimeSpacing,0]);
            end
        end
        
        function RXMatrix = GetRXMatrix(obj)
            % returns a matrix Q so that y=Q*r is the same as 
            % y=reshape(obj.Demodulation(r),[],1)     
            
            RXMatrix = zeros(obj.Nr.Subcarriers*obj.Nr.MCSymbols,obj.Nr.SamplesTotal);
            % could be implemented in a more efficient way, see "obj.GetSymbolNoisePower"
            for i_sample = 1:obj.Nr.SamplesTotal
                tempsignal = zeros(obj.Nr.SamplesTotal,1);
                tempsignal(i_sample) = 1;
                RXMatrix(:,i_sample) = reshape(obj.Demodulation(tempsignal),[],1);
            end
        end
        
        % Plot
        function [TransmitPower,Time] = PlotTransmitPower(obj, Rx)
            % plot the expected transmit power over time. The input 
            % argument represents the correlation of the data symbols. 
            % If no input argument is specified, an identity matrix is 
            % assumed (uncorrelated data)
            
            if exist('Rx','var')
                [V,D] = eig(Rx);
            else 
                % assume that Rx is an identity matrix, that is,
                % uncorrelated symbols
                V = eye(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
                D = V;               
            end
            D=sqrt(D);
            TransmitPower = zeros(obj.Nr.SamplesTotal,1);
            for i_lk = 1:obj.Nr.Subcarriers*obj.Nr.MCSymbols
                TransmitPower = TransmitPower+(abs(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2+abs(obj.Modulation(reshape(j*V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2)/2;
                if mod(i_lk,1000)==0
                    disp([int2str(i_lk/(obj.Nr.Subcarriers*obj.Nr.MCSymbols )*100) '%']);
                end
            end
            Time = (0:length(TransmitPower)-1)*obj.PHY.dt;
            if nargout==0
                plot(Time,TransmitPower);
                ylabel('Transmit Power');
                xlabel('Time(s)');
            end
        end
        
        function [PowerSpectralDensity,Frequency] = PlotPowerSpectralDensity(obj,Rx)
            % plot the power spectral density. The input argument 
            % represents the correlation of the data symbols. If no input
            % argument is specified, an identity matrix is assumed 
            % (uncorrelated data) 
            
            if exist('Rx','var')
                [V,D] = eig(Rx);
            else 
                V = eye(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
                D = V;               
            end
            D=sqrt(D);
            PowerSpectralDensity = zeros(obj.Nr.SamplesTotal,1);
            for i_lk = 1:obj.Nr.Subcarriers*obj.Nr.MCSymbols
                PowerSpectralDensity = PowerSpectralDensity+abs(fft(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk))).^2;
                if mod(i_lk,1000)==0
                    disp([int2str(i_lk/(obj.Nr.Subcarriers*obj.Nr.MCSymbols )*100) '%']);
                end
            end
            Frequency = (0:length(PowerSpectralDensity)-1)*1/(length(PowerSpectralDensity)*obj.PHY.dt);
            PowerSpectralDensity=PowerSpectralDensity/length(PowerSpectralDensity)^2/Frequency(2)^2;
            if nargout==0
                plot(Frequency,10*log10(PowerSpectralDensity));
                ylabel('Power Spectral Density (dB)');
                xlabel('Frequency (Hz)');
            end
        end
        
        function [PowerSpectralDensity,Frequency] = PlotPowerSpectralDensityUncorrelatedData(obj)
            % plot the power spectral density in case of uncorrelated data
            % symbols. Faster than FBMC.PlotPowerSpectralDensity
            
            PowerSpectralDensity = zeros(obj.Nr.SamplesTotal,1);
            for i_lk = 1:obj.Nr.Subcarriers
                V = zeros(obj.Nr.Subcarriers,obj.Nr.MCSymbols);
                V(i_lk,round(obj.Nr.MCSymbols/2))=1;
                PowerSpectralDensity = PowerSpectralDensity+abs(fft(obj.Modulation(V))).^2;
            end
            Frequency = (0:length(PowerSpectralDensity)-1)*1/(length(PowerSpectralDensity)*obj.PHY.dt);
            PowerSpectralDensity=obj.Nr.MCSymbols*PowerSpectralDensity/length(PowerSpectralDensity)^2/Frequency(2)^2;
            if nargout==0
                plot(Frequency,10*log10(PowerSpectralDensity));
                ylabel('Power Spectral Density (dB)');
                xlabel('Frequency (Hz)');
            end
        end
        
        % SIR and Noise power
        function Pn = GetSymbolNoisePower(obj, Pn_time)
            % returns the symbol noise power, that is, the noise power
            % after demodulation. The input argument is the noise power 
            % in the time domain. 
            
            RXFilterConvolutionMatrix = (convmtx(obj.Implementation.FilterImpulseResponseRX.',obj.Implementation.FFTSize));

            TempMatrix =  fft(conj(RXFilterConvolutionMatrix));
            TempMatrix2 = TempMatrix((1:obj.Nr.Subcarriers)+obj.Implementation.IntermediateFrequency,:);
            NoiseCorrlationMatrixNoEqualizer = sum(TempMatrix2.*conj(TempMatrix2),2);    % Samn as: diag(TempMatrix2*TempMatrix2')
            
            Pn = repmat(Pn_time*1/obj.Implementation.NormalizationFactor^2*NoiseCorrlationMatrixNoEqualizer.*abs(obj.Implementation.FilterEqualizerRX).^2,[1 obj.Nr.MCSymbols]);
            
        end
        
       function [SignalPower,InterferencePower] = GetSignalAndInterferencePowerQAM(...
                obj,...                                 % FOFDM object
                VectorizedChannelCorrelationMatrix,...  % Let the received signal be r=H*s with H representing a time-variant convolution matrix. Then "VectorizedChannelCorrelationMatrix" represents the expectation E{{H(:)*H(:)'}. We can obtain such matrix by ChannelModel.GetCorrelationMatrix
                DataSymbolCorrelationMatrix,...         % Correlation matrix of the vectorized data symbols
                TimeSamplesOffset,...                   % Time offset in samples (to improve the SIR)
                SubcarrierPosition,...                  % Subcarrier position for which the SIR is calculated.
                OFDMSymbolPosition...                   % OFDM symbol position in time for which the SIR is calculated.
                )
            % returns the signal and interference power for a
            % doubly-selective channel in case of QAM transmissions

            TXMatrix = obj.GetTXMatrix;
            RXMatrix = obj.GetRXMatrix;
            RXMatrix = [zeros(size(RXMatrix,1),TimeSamplesOffset),RXMatrix(:,1:end-TimeSamplesOffset)]; % Time offset compensation
              
            Index = SubcarrierPosition+(OFDMSymbolPosition-1)*obj.Nr.Subcarriers;

            % TempOld = kron(TXMatrix.',RXMatrix(Index,:))*VectorizedChannelCorrelationMatrix*kron(TXMatrix.',RXMatrix(Index,:))';
            % Much more efficient
            RXVectorRep = kron(sparse(eye(length(RXMatrix(Index,:)))),RXMatrix(Index,:)');
            Temp = RXVectorRep'*VectorizedChannelCorrelationMatrix*RXVectorRep;
            CorrMatrixNoData = TXMatrix.'*Temp*conj(TXMatrix);

            SignalDataSymbolCorrelationMatrix = zeros(size(DataSymbolCorrelationMatrix));
            SignalDataSymbolCorrelationMatrix(Index,Index) = DataSymbolCorrelationMatrix(Index,Index);
            InterferenceDataSymbolCorrelationMatrix = DataSymbolCorrelationMatrix;
            InterferenceDataSymbolCorrelationMatrix(Index,:)=0;
            InterferenceDataSymbolCorrelationMatrix(:,Index)=0;

            SignalPower = abs(sum(sum(CorrMatrixNoData.*SignalDataSymbolCorrelationMatrix)));
            InterferencePower = abs(sum(sum(CorrMatrixNoData.*InterferenceDataSymbolCorrelationMatrix)));
        end
       
        % Additional Functions
        function TimePos = GetTimeIndexMidPos(obj)
            % returns a vector which represents the discete time position
            % of the corresponding FBMC symbol (middle position)
            
            TimePos = obj.Implementation.ZeroGuardSamples+obj.Implementation.CyclicPrefix+round((obj.Implementation.CyclicPrefixLengthFilter+obj.Implementation.FilterLengthTX-1.5)/2)+round(obj.Implementation.FFTSize/2)+1+ (0:obj.Nr.MCSymbols-1)*obj.Implementation.TimeSpacing;
        end    
    end
end