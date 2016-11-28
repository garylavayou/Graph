classdef FlatClusterGraph < Graph
%% Cluster 1 is the high level cluster(core); the other clusters is the low level cluster.
    %   Detailed explanation goes here
    
    properties(Constant = true)
		NUM_PARALLEL_LINK = 3;
    end
    properties
		cluster_number;
		cluster_size;
	end
	
    methods
        function this = FlatClusterGraph()
        end
    end
    
    methods (Static)
        function this = Build(cluster_number, cluster_size, seed)
		% cluster_size is a scalar or vector
            if nargin>3
                rng(seed);
            else
                rng(12037);
            end
            this = FlatClusterGraph;
        
			if numel(cluster_size) == 1
				cluster_size = ones(cluster_number,1)*cluster_size;
			end
			this.cluster_size = cluster_size;
			this.cluster_number = cluster_number;
		    this.Size = sum(this.cluster_size);
            this.Adjacent = zeros(this.Size);

			% generate clusters
			for i = 1:cluster_number
				cluster_node_index = sum(this.cluster_size(1:i-1))+(1:this.cluster_size(i));
				this.Adjacent(cluster_node_index,cluster_node_index) = Graph.RandomAdjacent(this.cluster_size(i));
			end
			% connect clusters
			% generate virtual graph 0
			graph0 = Graph.RandomAdjacent(this.cluster_number);
			for i=1:this.cluster_number
				for j = 1:i-1
					if graph0(i,j) ~= 0
						l_n = min([NUM_PARALLEL_LINK, this.cluster_size(i), this.cluster_size(j)]);
						num_peer = randi(l_n);
						peer_index_i = sum(this.cluster_size(1:i-1))+unique_randi(this.cluster_size(i),num_peer, 'stable');
						peer_index_j = sum(this.cluster_size(1:j-1))+unique_randi(this.cluster_size(j),num_peer, 'stable');
						for k = 1:num_peer
							this.Adjacent(peer_index_i(k),peer_index_j(k)) = 1;
							this.Adjacent(peer_index_j(k),peer_index_i(k)) = 1;
						end
					end
				end
			end
            
            this.Renew;       % Update the properties in the object.
        end
    end
    
end

