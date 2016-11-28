%% FatTreeGraph
%   FatTreeGraph < Graph
%
% Data center topology.
% All the nodes have the same capacity, so all the links have the same capacity.
%%
classdef FatTreeGraph < Graph    
    properties
        N_core;
        N_pod;
        core_id;
        aggr_id;
        edge_id;
    end
    
    methods
        %% Constructor
        %   FatTreeGraph(K, C)
        %
        % |K| is a multiples of 2, |K|>=4.
        function this = FatTreeGraph(K, C)
            if nargin == 0 || isempty(K)
                return;
            end
            if nargin == 1
                C = 1;
            end
            total_size = (K/2)^2 + K^2;
            this.N_core = (K/2)^2;
            this.N_pod = K;
            this.Adjacent = zeros(total_size);
            this.core_id = 1:this.N_core;
            this.aggr_id = zeros(1, K^2/2);
            this.edge_id = zeros(1, K^2/2);
            %% Numbering and connecting
            % core nodes: 1:N_core
            % aggregation nodes in pod j=(1:K): N_core+(j-1)*N_pod+(1:K/2)
            % edge nodes in pod j=(1:K): N_core+(j-1)*N_pod+(K/2+1:K)
            for j = 1:this.N_pod
                this.aggr_id((j-1)*K/2+(1:K/2)) = this.N_core+(j-1)*this.N_pod+(1:K/2);
                this.edge_id((j-1)*K/2+(1:K/2)) = this.N_core+(j-1)*this.N_pod+(K/2+1:K);
                for k = 1:K/2
                    %%%
                    % the k'th node in the aggregation layer connect to the core nodes with
                    % number (k-1)*K/2+(1:K/2);
                    for i = (k-1)*K/2+(1:K/2)
                        this.Adjacent(this.N_core+(j-1)*this.N_pod+k,i) = C;
                    end
                    for i = this.N_core+(j-1)*this.N_pod+(K/2+1:K)
                        this.Adjacent(this.N_core+(j-1)*this.N_pod+k,i) = C;
                    end
                end
            end
            
            this.Adjacent = this.Adjacent + this.Adjacent';
            if ~Graph.VerifyConnectivity(this.Adjacent)
                error('Error: Graph is not all-connected.');
            end
            
            this.Renew;      % Update the properties in the object.
            
        end

        function Copy(this, graph)
            Copy@Graph(this, graph);        % call the method of superclass
            this.N_core = graph.N_core;
            this.N_pod = graph.N_pod;
            this.core_id = graph.core_id;
            this.aggr_id = graph.aggr_id;
            this.edge_id = graph.edge_id;
        end
    end

end

