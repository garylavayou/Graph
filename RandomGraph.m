classdef RandomGraph < Graph
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function this = RandomGraph(N, seed)
            if nargin == 0 || isempty(N)
                return;
            end            
            if nargin>1
                rng(seed);
            else
                rng(12037);
            end
            
            this.Adjacent = RandomGraph.RandomAdjacent(N);
            
            this.Renew;           % Update the properties in the object.     
        end
    end
    
end

