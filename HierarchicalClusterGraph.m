classdef HierarchicalClusterGraph < Graph
    %HIERARCHICALCLUSTERGRAPH 2-level clustered graph
    %   Cluster 1 is the high level cluster(core) and the other clusters are the low level
    % cluster (RAN cluster). 
    
    properties(Constant = true)
		NUM_CANDIDATE_CORE_NODE = 2;
    end
    properties
		cluster_number;
		cluster_size;
	end
	
    methods
        function this = HierarchicalClusterGraph()
        end
        
        function SetLinkCapacity(this, links, capacity)
            % links: index of links. If links is a vector, it represent the numbering of
            %     links. if links is a n*2 matrix, it represents the (head tail) form of
            %     links. if links is a length 2 cell array , it represents the subgraph's
            %     node numbering.
            % capacity: if capacity is a scalar, all links have the capacity; if capacity
            % is a vector, it must have the length that equal to the number of links.
            % 
            if iscell(links)
                if length(links) ~= 2
                    error('error: Length of cell array LINKS is not 2.');
                end
                if isscalar(capacity)
                    this.Adjacent(links{1}, links{2}) = ...
                        this.Adjacent(links{1}, links{2}) * capacity;
                else
                    %% TODO
                end
            elseif isvector(links)
                %% TODO
            else 
                if size(links,2) == 2  % link (i,j);
                    %% TODO
                elseif size(links,1) == 2 % links in (i1:i2, j1:j2)

                else
                    error('error: the size of input argument links is invalid.');
                end
            end
            % update the Capacity
            [~, ~, this.Capacity] = find(this.Adjacent);
        end
    end
    
    methods (Static)
        function this = Build(cluster_number, cluster_size, seed, link_capacity)
            % cluster_number: Number of clusters.
            % cluster_size: Size of a cluster. cluster_size is a scalar or vector. When
            %     cluster is a scalar, all clusters have the same size.
            % link_capacity: specifies the link capacity of the core cluster, the RAN
            %     cluster, and the link between core and RAN clusters. 
            %         if link_capacity is a scalar, all links have the same capacity.
            %         if link_capacity is a vector with length 3, then core cluster, the
            %         RAN cluster, and the link between core has different capacity
            %         parameter. 
            %         if link_capacity is a vector with length (cluster_number+1), then
            %         each RAN cluster has different capacity parameter.
            %     NOTE: the link capacity with one RAN cluster can be different according
            %         its distance to the D-GW. We do not consider it at present.

            if nargin>=3 && isempty(seed) == 0
                rng(seed);
            else
                rng(12037);
            end
            if nargin>=4 && isempty(link_capacity)==0
                if isscalar(link_capacity)
                    link_capacity = ones(cluster_number+1,1) * link_capacity;
                elseif isvector(link_capacity)
                    if length(link_capacity) == 3
                        link_capacity = [link_capacity(1); ...
                            link_capacity(2)*ones(cluster_number-1,1); link_capacity(3)];
                    elseif length(link_capacity) ~= cluster_number+1
                        error('error: the length of input argument link_capacity is invalid.');
                    else
                        %% TODO
                    end
                else
                    error('error: the size of input argument link_capacity is invalid.');
                end
            else
            end
            
            this = HierarchicalClusterGraph;
        
			if numel(cluster_size) == 1
				cluster_size = ones(cluster_number,1)*cluster_size;
			end
			this.cluster_size = cluster_size;
			this.cluster_number = cluster_number;
            total_size = sum(this.cluster_size);
            this.Adjacent = zeros(total_size);

			% generate clusters
			for i = 1:cluster_number
				node_index = sum(this.cluster_size(1:i-1))+(1:this.cluster_size(i));
				this.Adjacent(node_index,node_index) = Graph.RandomAdjacent(this.cluster_size(i));
                if nargin >= 4
                    if isempty(link_capacity)
                        %% TODO: hierarchical capacity in the RAN as SD-RAN
                    else
                        this.Adjacent(node_index,node_index) = ...
                            this.Adjacent(node_index,node_index)*link_capacity(i);
                    end
                else
                    %% TODO: link capacity is proportion to the network size (this.graph.Size)
%                     if i == 1
%                         this.Adjacent(node_index,node_index) = ...
%                             this.Adjacent(node_index,node_index)*100;
%                     end
                end
			end
			% connect core and clusters
			core_size = this.cluster_size(1);
			core_node_index = 1:core_size;
			core_graph = this.Adjacent(core_node_index,core_node_index);
            for i = 2:cluster_number
                default_core_node_id = randi(core_size);
                core_node_neighbor = find(core_graph(default_core_node_id,:)~=0);
                l_n = min(HierarchicalClusterGraph.NUM_CANDIDATE_CORE_NODE, length(core_node_neighbor));
                num_peer = randi([0 l_n])+1;
                core_peer_index = [default_core_node_id ...
                    core_node_neighbor(unique_randi(length(core_node_neighbor), num_peer-1, 'stable'))];
                cluster_peer_index = sum(this.cluster_size(1:i-1))+...
                    unique_randi(this.cluster_size(i), num_peer, 'stable');
                for j = 1:num_peer
                    % TODO: link capacity is proportion to the network size (this.graph.Size)
                    
                    if nargin >= 4
                        if isempty(link_capacity)
                            %% TODO: hierarchical capacity in the RAN as SD-RAN
                        else
                            this.Adjacent(core_peer_index(j),cluster_peer_index(j)) = ...
                                link_capacity(cluster_number+1);
                            this.Adjacent(cluster_peer_index(j),core_peer_index(j)) = ...
                                link_capacity(cluster_number+1);
                        end
                    else
                        %% TODO: link capacity is proportion to the network size (this.graph.Size)
                        this.Adjacent(core_peer_index(j),cluster_peer_index(j)) = 1;
                        this.Adjacent(cluster_peer_index(j),core_peer_index(j)) = 1;
                    end
                end
            end
            
            if Graph.VerifyConnectivity(this.Adjacent) == 0
                error('error: The graph is not all connected.');
            end
            
            this.Renew;       % Update the properties in the object.

        end
    end
    
end

