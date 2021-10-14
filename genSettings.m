classdef genSettings < handle
    %this is a metaclass that defines global settings for all other classes
    
    properties (Constant)
        %p-p-1 convention is the one used in Plethora (the corner belongs
        %to K_p union bar K_{p-1}
        %p-p convention means the corner belongs to K_p union bar K_p
        %in the software, only xill and xirr functions are available in p-p
        %convention; I thought it's simpler, but it's actually equivalent
        %in difficulty
        
        %CALC_CONVENTION = 'p-p'; 
        CALC_CONVENTION = 'p-p-1';
        
        CALC_SHOW_MESSAGE = false; %shows which terms are non-zero in calcL
        
        EMPTY_BLOCKS = false; %Plethora defines additionally intervals for empty blocks
                              %false means we don't define them
                              
    end
    
    methods
        function obj = settings()
            
        end
        
    end
end

