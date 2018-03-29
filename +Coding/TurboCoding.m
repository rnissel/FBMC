classdef TurboCoding < handle 
  % =====================================================================        
  % This MATLAB class represents a turbo coder (LTE).
  % It requires the MATLAB Communications System Toolbox!
  % Usage: 
  %    1) CodingObject = Coding.TurboCoding(NrDataBits,NrCodedBits)
  %    2) CodingObject.TurboEncoder(Bitstream)
  %    3) CodingObject.TurboDecoder(LLR)
  % Additionally, we might want to update the internal interleaver by
  %    CodingObject.UpdateInterleaving; 
  % The code is based on the book 
  % "Understanding LTE with MATLAB", Houman Zarrinkoub
  % =====================================================================    
  % Ronald Nissel, rnissel@nt.tuwien.ac.at
  % (c) 2017 by Institute of Telecommunications, TU Wien
  % www.nt.tuwien.ac.at
  % =====================================================================    

  properties (SetAccess = private)
      CommTurboEncoder
      CommTurboDecoder
      NrDataBits
      NrCodedBits
      CodeRate
      Interleaving
  end
  
  methods
  	function obj = TurboCoding(varargin)
        if varargin{2}>varargin{1}
            obj.NrDataBits = varargin{1}; 
            obj.NrCodedBits = varargin{2};
        else
            obj.NrDataBits = varargin{2}; 
            obj.NrCodedBits = varargin{1};
        end
        obj.CodeRate = obj.NrDataBits/obj.NrDataBits;
        
        obj.Interleaving = randperm(obj.NrDataBits);

        obj.CommTurboEncoder = comm.TurboEncoder('TrellisStructure', poly2trellis(4, [13 15], 13),'InterleaverIndicesSource','Input port');    
        obj.CommTurboDecoder = comm.TurboDecoder('TrellisStructure', poly2trellis(4, [13 15], 13),'InterleaverIndicesSource','Input port', 'NumIterations', 10);
       
    end  
        
    function CodedBits = TurboEncoder( obj , DataBits )
        ImplementationCodeRate = (obj.NrDataBits+4)/obj.NrCodedBits;

        CodedBits_1_3 = step(obj.CommTurboEncoder,DataBits,obj.Interleaving);

        if ImplementationCodeRate>=1/3
            CodedBits=RateMatcher(CodedBits_1_3,length(DataBits),ImplementationCodeRate);
        else
            IndexTurboRepetition = repmat(1:length(CodedBits_1_3),1,ceil(obj.NrCodedBits/length(CodedBits_1_3)));
            IndexTurboRepetition=IndexTurboRepetition(1:obj.NrCodedBits);
            CodedBits = CodedBits_1_3(IndexTurboRepetition);
        end
    end   
    
    function UpdateInterleaving(varargin)
        obj = varargin{1};
        if numel(varargin)>=2            
            obj.Interleaving = varargin{2};              
        else
            obj.Interleaving = randperm(obj.NrDataBits);
        end
    end
    
    function DecodedDataBits = TurboDecoder( obj , LLR ) 
        NrCodedBits_1_3 = (obj.NrDataBits+4)*3;
        ImplementationCodeRate = (obj.NrDataBits+4)/obj.NrCodedBits;


        if ImplementationCodeRate>=1/3
            LLR_RateMatched = RateDematcher(LLR, obj.NrDataBits);
        else
            IndexTurboRepetition = repmat(1:NrCodedBits_1_3,1,ceil(obj.NrCodedBits/NrCodedBits_1_3));
            IndexTurboRepetition=IndexTurboRepetition(1:obj.NrCodedBits);

            LLR_RateMatched = Coding.LTE_common_turbo_rate_matching_bit_selection_and_pruning(LLR,...
                IndexTurboRepetition-1,...
                2,...
                NrCodedBits_1_3).';
        end
        DecodedDataBits = step(obj.CommTurboDecoder,LLR_RateMatched,obj.Interleaving);
        end
    
  end
end


function CodedBits = RateMatcher(CodedBits_1_3, NrDataBits, CodeRate)
    NrDataBitsP4 = NrDataBits+4;
    a = 32;
    b = ceil(NrDataBitsP4/a);

    PermutationPattern      = [0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30, 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31];
    Temp                    = reshape([nan(1,a*b-NrDataBitsP4),1:NrDataBitsP4],a,b).';
    Temp2                   = reshape(Temp.',[],1);
    Index1   = reshape(Temp(:,PermutationPattern+1),[],1);
    Index2   = Temp2(PermutationPattern(floor((0:a*b-1)/b)+1)+a*mod(0:a*b-1,b)+1);

    c0      = SubBlockInterleaving(CodedBits_1_3(1:3:end),Index1);
    c12     = reshape([SubBlockInterleaving(CodedBits_1_3(2:3:end),Index1),SubBlockInterleaving(CodedBits_1_3(3:3:end),Index2)].',[],1);

    c = [c0(isfinite(c0));c12(isfinite(c12))];

    CodedBits = c(1:round(NrDataBitsP4/CodeRate));
end

function LLR_RateMatched = RateDematcher(LLR, NrDataBits)
    NrDataBitsP4 = NrDataBits+4;
    a = 32;
    b = ceil(NrDataBitsP4/a);

    PermutationPattern      = [0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30, 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31];
    Temp                    = reshape([nan(1,a*b-NrDataBitsP4),1:NrDataBitsP4],a,b).';
    Temp2                   = reshape(Temp.',[],1);
    Index1   = reshape(Temp(:,PermutationPattern+1),[],1);
    Index2   = Temp2(PermutationPattern(floor((0:a*b-1)/b)+1)+a*mod(0:a*b-1,b)+1);

    LLR0 = zeros(3*NrDataBitsP4,1);
    LLR0(1:numel(LLR)) = LLR;  
    
    l0 = SubBlockDeInterleaving(LLR0(1:NrDataBitsP4), Index1);
    l12 = reshape(SubBlockDeInterleaving(LLR0(NrDataBitsP4+1:end), reshape([Index1,Index2+NrDataBitsP4].',[],1)), NrDataBitsP4 , 2);
    
    LLR_RateMatched = reshape([l0 l12].',[],1);
end

function c = SubBlockInterleaving(c_in,Index)
    c                   = zeros(size(Index));
    c(~isnan(Index)==1) = c_in(Index(~isnan(Index)==1));
    c(isnan(Index)==1)  = nan*ones(sum(isnan(Index)==1),1);
end


function l = SubBlockDeInterleaving(LLR0,Index)
    l                           = zeros(size(LLR0));
    l(Index(~isnan(Index)==1))  = LLR0;
end
  
      