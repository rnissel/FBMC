% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2017 by Institute of Telecommunications, TU Wien
% www.tc.tuwien.ac.at    

% This function calculates the bit error probability for an arbitrary
% signal constellations in an AWGN channel
function BitErrorProbability = BitErrorProbabilityAWGN(...
    SNR_dB, ...       % The SNR in the complex domain => SNR_FBMC = SNR_OFDM-3dB.
    SymbolMapping, ...% The symbol mapping with mean(SymbolMapping.*conj(SymbolMapping))==1. For example in 4QAM we have: SymbolMapping=[0.7071 + 0.7071i;-0.7071 + 0.7071i;0.7071 - 0.7071i;-0.7071 - 0.7071i];
    BitMapping)       % The bitmapping corresponding to the symbol mapping. For example for 4QAM we have: BitMapping = [0 0;1 0;0 1;1 1];



    % For the decision regions we assume a rectangular regular grid. The rest
    % of the function could also be used for an irregular grid but
    % the decision regions would have to be rewritten!        
    HalfedDecisionInterval = min(abs(real(SymbolMapping)));
    DecisionRegions = [...
        real(SymbolMapping)- HalfedDecisionInterval ...
        real(SymbolMapping)+ HalfedDecisionInterval ...
        imag(SymbolMapping)- HalfedDecisionInterval ...
        imag(SymbolMapping)+ HalfedDecisionInterval  ];
    DecisionRegions(min(real(SymbolMapping))==real(SymbolMapping),1) = -inf;
    DecisionRegions(max(real(SymbolMapping))==real(SymbolMapping),2) = +inf;
    DecisionRegions(min(imag(SymbolMapping))==imag(SymbolMapping),3) = -inf;
    DecisionRegions(max(imag(SymbolMapping))==imag(SymbolMapping),4) = +inf;
    
    BitErrorProbability = nan(length(SNR_dB),1);
    for i_SNR = 1:length(SNR_dB)
        Pn = 10^(-SNR_dB(i_SNR)/10);
        ProbabilityMatrix = nan(size(SymbolMapping,1),size(SymbolMapping,1));
        for i_symbol = 1:size(SymbolMapping,1)
            x = SymbolMapping(i_symbol);
            ProbabilityMatrix(:,i_symbol)=GaussianRatioProbabilityRectangularRegion(x,Pn,DecisionRegions(:,1),DecisionRegions(:,2),DecisionRegions(:,3),DecisionRegions(:,4));
        end 
        ErrorProbability = nan(2,size(BitMapping,2));
        for i_bit= 1:size(BitMapping,2)
            for i_zero_one = [0 1]
                index_x = (BitMapping(:,i_bit)==i_zero_one);
                ErrorProbability(i_zero_one+1,i_bit) = mean(sum(ProbabilityMatrix(not(index_x),index_x)));
            end   
        end
        BitErrorProbability(i_SNR) = mean(mean(ErrorProbability));
    end
end


% This function calculates the Probability that y=x+n is within the 
% rectrangular region "(zRlower zIlower] \times (zRupper zIupper]". 
% It requires the function "GaussianRatioCDF"
function Probability = GaussianRatioProbabilityRectangularRegion(...
    x,...               % Transmitted symbol
    Pn,...              % Complex noise power
    zRlower,...         % Determines the rectangular region
    zRupper,...
    zIlower,...
    zIupper)


    CDF_RegionA = GaussianCDF(x,Pn,zRupper,zIupper);
    CDF_RegionB = GaussianCDF(x,Pn,zRlower,zIlower);
    CDF_RegionC = GaussianCDF(x,Pn,zRlower,zIupper);
    CDF_RegionD = GaussianCDF(x,Pn,zRupper,zIlower);
    
    Probability = CDF_RegionA+CDF_RegionB-CDF_RegionC-CDF_RegionD;
    
end



% This function calculates the CDF, that is Pr(real(x)<zR & imag(x)<zR), 
% whereas x is a complex Gaussian random variable
function CDF = GaussianCDF(...
    x,...               % ....
    Pn,...              % ...
    zR,...              % Real part of the CDF, i.e, Pr(real(y/h)<zR &...)
    zI)                 % Imaginary part of the CDF, i.e, Pr(... & imag(y/h)<zR)


    CDF = normcdf(zR,real(x),sqrt(Pn/2)).*normcdf(zI,imag(x),sqrt(Pn/2));

end