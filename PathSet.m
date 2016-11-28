classdef PathSet < handle
    %PATHSET
    
    properties(Constant)
        L0 = 4;
    end
    
    properties
        % PATH_SET (N*N) cell array, in which each element {i,j} is also a cell array,
        %     which reprensent the list of path from i to j. 
        Element;
        % NUMBER number of path in the path set.
        Number;
        MaxId;
        
        % REVERSE_INDEX
        reverse_index;
        max_id;
    end
    methods
        function n = get.Number(this)
            n = this.Number;
        end
        function n =get.MaxId(this)
            n = this.max_id;
        end
    end
    
    methods
        function this = PathSet(n1,n2)
			switch nargin
				case 1
					this.Element = cell(n1);
				case 2
					this.Element = cell(n1,n2);
				otherwise
					error('error: invalid number of arguments.');
			end
            this.Number = 0;
        end
        function b = IsEmpty(this)
            b = this.Number==0;
        end
		
        function new_path_index = Add(this, path_list, bandwidth)
            %ADD add PATH_LIST into the PATH_SET.
            % path_list: the list of path to be added.
			%    path_list is a cell array whose element is a vector of nodes that
			%    comprise a path. 
			% bandwidth: the bandwidth of paths in the path list. 
			%     if bandwidth is a scalar, all path have the same bandwidth
			%     if bandwidth is a vector, it has the same length with path_list.
			%     if bandwidth is not specified, the paths have default bandwidth 1.
			if nargin >=3 && isempty(bandwidth)==0
				if isscalar(bandwidth)
					bandwidth = bandwidth * ones(length(path_list),1);
				elseif isvector(bandwidth)
					if length(bandwidth) ~= length(path_list)
						error('error: Input arguments BANDWIDTH and PATH_LIST should have the same length.');
					end
				else
					error('error: Input argument BANDWIDTH must be a vector or scalar, otherwise do not set this argument.');
				end
			else
				bandwidth = ones(length(path_list),1);
			end
            new_path_index = zeros(length(path_list),3);
            for i = 1:length(path_list)
                path = Path(path_list{i}, bandwidth(i));
                src = path.Head;
                dest = path.Tail;
                index = Index(this, path, src, dest);
                if isempty(index)   % no path exists
                    this.Element{src,dest} = path;
                    new_path_index(i,:)=[src dest 1];
                    this.Number = this.Number + 1;
                elseif sum(index)==0  % this path is different from existing path
                    path_num = length(this.Element{src,dest});
                    this.Element{src,dest}(path_num+1)=path;
                    new_path_index(i,:)=[src, dest, path_num+1];
                    this.Number = this.Number + 1;
                else % this path already exist
					% increase the bandwidth
                    new_path_index(i,:)=index;
					this.Element{src,dest}(index(3)).bandwidth = ...
						this.Element{src,dest}(index(3)).bandwidth + path.bandwidth;
                end
            end
        end
		
        function new_path_index = AddPath(this, path, src_index, dest_index)
            %ADDPATH add PATH_LIST into the PATH_SET.
            % path: the path to be added.
            % src_index: source node's index in the path_set.
            % dest_index: destination node's index in the path_set.
            % NOTE: the src_index and dest_index is essential, since in some network(RadioAccessNetwork) the path_set may not use the node id as index so as to compress storage. In other words, if the path_set is indexed by node id, these two arguments can be elide.
            if nargin == 2
                src_index = path.Head;
                dest_index = path.Tail;
            end
            index = Index(this, path, src_index, dest_index);
            if isempty(index)   % no path exists
                this.Element{src_index,dest_index} = Path(path.node_list, path.bandwidth);
                new_path_index=[src_index dest_index 1];
                this.Number = this.Number + 1;
            elseif sum(index)==0  % this path is different from existing path
                path_num = length(this.Element{src_index,dest_index});
                this.Element{src_index,dest_index}(path_num+1)=...
                    Path(path.node_list, path.bandwidth);
                new_path_index=[src_index dest_index path_num+1];
                this.Number = this.Number + 1;
            else % this path already exist
                new_path_index=index;
            end
        end
        function index = Index(this, path, src_index, dest_index)
            %PATH_INDEX retrive the index of a PATH in the PATH_SET.
            % path: request path reprensented by an vector.
            if IsEmpty(this)
                index = [];
            else
                path_num = length(this.Element{src_index,dest_index});
                b_equal = 0;
                for j=1:path_num
                    if path.EqualPath(this.Element{src_index,dest_index}(j))
                        b_equal = 1;
                        break;
                    end
                end
                if ~b_equal
                    index = [0 0 0];
                else
                    index = [src_index,dest_index,j];
                end
            end
        end
        function path_num = Count(this, src_index, dest_index)
            if IsEmpty(this)
                path_num = 0;
            else
                path_num = length(this.Element{src_index,dest_index});
            end
        end
        function AllocateId(this, L)
            % path is indexed by column
            % if L is given, only allocate id to former L paths.
            id = 0;
            this.reverse_index = zeros(this.Number, 3);
            for j = 1:size(this.Element,2)
                for i = 1:size(this.Element,1)
                    if ~isempty(this.Element{i,j})
                        if nargin == 2
                            num_path = min(L, length(this.Element{i,j}));
                        else
                            num_path = length(this.Element{i,j});
                        end
                        for k=1:num_path
                            id = id + 1;
                            this.Element{i,j}(k).id = id;
                            this.reverse_index(id,:) = [i j k];
                        end
                    end
                end
            end
            this.max_id = id;
        end
        
        function Copy(this, ps)
			[s1,s2] = size(ps.Element);
			this.Element = cell(s1,s2);
			for i = 1:s1
				for j = 1:s2
					n = length(ps.Element{i,j});
					for k = 1:n
						path = ps.Element{i,j}(k);
						this.Element{i,j}(k) = Path(path.node_list, path.bandwidth);
					end
				end
			end
				
            this.reverse_index = ps.reverse_index;
            this.max_id = ps.max_id;
        end
		
		function f = MultiplyBandwidth(this, factor)
            % f: f(i,j) represents total flow amount from node i to node j;
			[s1,s2] = size(this.Element);
			f = zeros(s1,s2);
			for i = 1:s1
				for j = 1:s2
					n = length(this.Element{i,j});
					if n ~= 0
						for t = 1:n
							this.Element{i,j}(t).MultiplyBandwidth(factor);
							f(i,j) = f(i,j) +  this.Element{i,j}(t).bandwidth;
						end
					end
				end
			end
		end
		
		function fe = CountEdgeFlow(this, graph, factor)
			if nargin<=2
				factor = 1;
			end
			fe = zeros(length(graph.Capacity),1);
            [s1,s2] = size(this.Element);
			for i = 1:s1
				for j = 1:s2
					n = length(this.Element{i,j});
					if n == 0
						continue;
					end
					for t = 1:n
						if factor~=1
							this.Element{i,j}(t).MultiplyBandwidth(factor);
						end
						p = this.Element{i,j}(t).node_list;    % =>dest_set(i) => flow(i)
						for r=1:length(p)-1;
							ei = graph.Inverse(p(r),p(r+1));
							fe(ei) = fe(ei) + this.Element{i,j}(t).bandwidth;
						end 
					end
				end
			end
		end
    end
    methods (Static)
        % TODO move this function to the radio_access_network class.
        function path_set = Build(network)
            % Note: flow_list & user_bs corresponds to graph
            graph = network.ran_graph;
            user_bs = network.user_bs(network.active_user_index);
            router_id = network.index.Router;
            num_wire_nodes = graph.NumberWireNode;
            
            Gl = full(graph.Adjacent(1:num_wire_nodes, 1:num_wire_nodes));
            nz = find(Gl~=0);
            Gl(nz) = 1./Gl(nz);
            Gl(Gl==0) = inf;
            path_set = PathSet(graph.NumberRouter, graph.NumberUser);
            % find multiple path from router to user.
            for i = 1:graph.NumberRouter
                path_list = Graph.DijkstraTree(Gl, router_id(i), 1:graph.NumberBS);
                user_list = network.flow_list(network.flow_list(:,1) == router_id(i), 2);
                for j = 1:length(user_list)
                    % when indexing fake_user_bs, the index should be 1-based.
                    user_index = user_list(j) - network.UserIdOffset;
                    [~,ix] = sort(user_bs{user_index}(:,2));
                    candidate_BS = user_bs{user_index}(ix(1:min(PathSet.L0,length(ix))),1);
                    for k = 1:length(candidate_BS)
                        % denoting the graph elements, user index is the user id in the
                        % flow table.
                        path = Path([path_list{candidate_BS(k)} user_list(j)], 0);
                        path_set.AddPath(path, i, user_index);
                    end
                end
            end
            path_set.AllocateId();
        end
    end
end

