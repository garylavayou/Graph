classdef CoreBackhaulGraph < Graph
    %COREBACKHAULGRAPH 2-layer graph comprises core node and backhaul aggregate node
    %   the RAN cluster is abstracted as an aggregate node.
    
    properties
        AccessCoreMap;
        NumCoreNodes;
        NumCoreEdges;
    end
    
    methods (Access=private)
        function this = CoreBackhaulGraph()
        end
    end
    
    methods (Static)
        %% Function: build graph
        function this = Build(N_core, var_an_num, cdf_an_num, var_an_deg, cdf_an_deg, seed)
            %% Random Graph
            % var_an_num: the range of number of access node connect to one core node
            % cdf_an_num: corresponding cdf of var_an_num
            var_an_deg = var_an_deg-1;
            if nargin > 5
                rng(seed);
            else
                rng(36590);
            end
            this = CoreBackhaulGraph;
            
            G_core = Graph.RandomAdjacent(N_core);
            
            %% specify the association between core nodes and access nodes.
            access_nodes = cell(1, N_core);
            max_an_num = var_an_num(end);   % max(var_an_num);
            id = N_core;
            for i =1:N_core
                % determine the number of access nodes connected to core node i.
                an = zeros(1,max_an_num);
                l_an = CoreBackhaulGraph.rand_cdf(var_an_num, cdf_an_num);
                % numbering the access nodes
                for j = 1:l_an
                    id =id + 1;
                    an(j) = id;
                end
                access_nodes{i} = an(an~=0);
            end
            this.Size = id;
            this.Adjacent = zeros(id);
            this.Adjacent(1:N_core,1:N_core) = G_core;
            access_core_map = zeros(this.Size - N_core, 1);
            for i = 1:N_core
                for j = 1:length(access_nodes{i})
                    id = access_nodes{i}(j);
                    this.Adjacent(i,id) = inf;   % default association, backhaul link has infinite capacity.
                    access_core_map(id-N_core) = i;
                    neighbor_i = find(G_core(i,:)~=0);
                    l_n = min(CoreBackhaulGraph.rand_cdf(var_an_deg, cdf_an_deg), length(neighbor_i)); % number of neighbours
                    if l_n == 0
                        continue
                    end
                    new_index = unique(randi(length(neighbor_i),2*l_n,1),'stable');
                    new_index = new_index(1:min(l_n,length(new_index)));
                    for k=new_index'
                        this.Adjacent(neighbor_i(k),id) = inf;
                    end
                end
            end
            fprintf('2-layer graph information:\n');
            fprintf('\t average degree: %f\n',mean(sum(this.Adjacent~=0)));
            %% limit the backhaul link's capacity
            % By default, we can set the backhaul link capacity to INFINITY, thus giving the
            % full flexibility to transfer traffic.
            % Alternatively, we can limit the backhaul link capacity to simulate the real
            % network scenario. Considering the core node, its capacity in the core network
            % side should be balanced with the access network side. hence, we let the access
            % link fairly share the total capacity to the core network. Additionally, core
            % nodes may have traffic as well, such as content distribution.
            %         core_degree = (sum(G_core,1))';            % in
            %         core_access_degree = sum(Gc(1:N_core,N_core+1:N)~=0,2);  % out
            %         access_capacity = core_degree./(core_access_degree+1);
            % Since INF in MATLAB is difficult to handle, we deal with the default case by
            % setting the capacity with a relatively large value. there, we set it to
            % number of nodes (N).
            access_capacity = this.Size*ones(N_core,1);
            for i = 1:N_core
                this.Adjacent(i,this.Adjacent(i,:)==inf)=access_capacity(i);
            end
            this.AccessCoreMap = access_core_map;		% mapping the default association.
            this.NumCoreNodes = N_core;
            this.NumCoreEdges = length(find(G_core~=0));
            
            this.Renew;      % Update the properties in the object.   
        end
            
        %% generate random number according to the distribution.
        function p = rand_cdf(var, cdf)
            p = var(end);
            x =rand;
            for i = 1:length(var)
                if x < cdf(i)
                    p = var(i);
                    break;
                end
            end
        end
        
    end
    
end

