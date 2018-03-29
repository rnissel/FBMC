classdef SpaceCoding < handle 
    % =====================================================================        
    % This MATLAB class represents a basic implementation of Alamouti's 
    % space time(or frequency) block code.
    % =====================================================================    
    % Ronald Nissel, rnissel@nt.tuwien.ac.at
    % (c) 2017 by Institute of Telecommunications, TU Wien
    % www.tc.tuwien.ac.at
    % =====================================================================    

    properties (SetAccess = private)
        Method
        RowSpreading
        Nr
    end
  
    methods
        function obj = SpaceCoding(Input_Method,Input_RowSpreading)
            obj.Method          = Input_Method;  
            obj.RowSpreading    = Input_RowSpreading; 
            switch obj.Method
                case 'Alamouti2x1'
                    obj.Nr.TransmitAntennas = 2;
                    obj.Nr.ReceiveAntennas  = 1;                
                otherwise
                    error('SpaceCode not supported!')
        	end              
        end  
        
        function CodedDataSymbols = Encoder(obj,DataSymbols)
        
            if obj.RowSpreading
                DataSymbols = permute(DataSymbols,[2 1 3]);
            end
        
            switch obj.Method
                case 'Alamouti2x1'
                    CodedDataSymbols                = nan([size(DataSymbols),obj.Nr.TransmitAntennas]);
                    CodedDataSymbols(1:2:end,:,1)   = DataSymbols(1:2:end,:);
                    CodedDataSymbols(2:2:end,:,1)   = -conj(DataSymbols(2:2:end,:));
                    CodedDataSymbols(1:2:end,:,2)   = DataSymbols(2:2:end,:);
                    CodedDataSymbols(2:2:end,:,2)   = conj(DataSymbols(1:2:end,:));                             
                otherwise
                    error('SpaceCode not supported!')
            end  

            if obj.RowSpreading
                CodedDataSymbols = permute(CodedDataSymbols,[2 1 3]);
            end

        end 

        function DecodedDataSymbols = Decoder(obj,ReceivedDataSymbols,Channel)

            if obj.RowSpreading
                ReceivedDataSymbols = permute(ReceivedDataSymbols,[2 1 3]);
                Channel             = permute(Channel,[2 1 3]);
            end

            switch obj.Method
                case 'Alamouti2x1'
                    ReceivedDataSymbols(2:2:end,:) = conj(ReceivedDataSymbols(2:2:end,:));

                    ReceivedDataSymbolsShaped   = reshape(ReceivedDataSymbols,2,[]);
                    Channel1Reshaped            = mean(reshape(Channel(:,:,1),2,[]),1);
                    Channel2Reshaped            = mean(reshape(Channel(:,:,2),2,[]),1);

                    DecodedDataSymbols  = nan(size(ReceivedDataSymbolsShaped));
                    for i_index = 1:size(ReceivedDataSymbolsShaped,2)
                        DecodedDataSymbols(:,i_index) = 2/((abs(Channel1Reshaped(i_index))^2+abs(Channel2Reshaped(i_index))^2)*sqrt(2))*...
                        [ Channel1Reshaped(i_index),Channel2Reshaped(i_index);...
                          conj(Channel2Reshaped(i_index)),-conj(Channel1Reshaped(i_index)) ]'...
                          * ReceivedDataSymbolsShaped(:,i_index); 
                    end
                    DecodedDataSymbols = reshape(DecodedDataSymbols,size(ReceivedDataSymbols));            
                otherwise
                    error('SpaceCode not supported!')
            end   

            if obj.RowSpreading
                DecodedDataSymbols = permute(DecodedDataSymbols,[2 1 3]);
            end
        end       


    end
end

      