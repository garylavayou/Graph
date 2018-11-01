%% Directed Graph
% directed graph
%% TODO
% Refactor the adjacent matrix: rows as tails and columns as head of edges, so
% as to keep consistent with graph/digraph class.

%%
classdef DirectedGraph < matlab.mixin.Copyable
	properties (SetAccess=protected)
		Adjacent;       % matrix link weight (TODO: differentiate from Delay)
		Capacity;       % matrix link capacity
		
		Head;
		Tail;
		LinkWeight;
		%         LinkCapacity;
	end
	
	properties (Access = private)
		link_id;         % inverse index of links
		mat_cost;        % for APSP
		mat_path;        % for APSP
		mat_cpath;       % for APSP
	end
	properties(Dependent)
		% EdgeTable;
		NumberNodes;
		NumberEdges;
	end
	
	methods
		%% Constructor
		function this = DirectedGraph(A, C)
			switch nargin
				case 1
					if isa(A, 'DirectedGraph')
						this = A.copy;
					elseif isa(A, 'digraph')
						constructByEdgeNodeTable(A.Edges, A.Nodes);
					elseif isnumeric(A)
						constructByAdjacent(A, []);
					elseif isstruct(A) 
						if isfield(A, 'Edges') && istable(A.Edges)
							% A includes an EdgeTable (and a NodeTable).
							edge_table = A.Edges;
							if isfield(A, 'Nodes')
								node_table = A.Nodes;
							else
								node_table = [];
							end
							constructByEdgeNodeTable(edge_table, node_table);
						elseif isfield(A, 'Adjacent')
							% A includes an adjacency matrix (and the corresponding capacity matrix).
							Adj = A.Adjacent;
							if isfield(A, 'Capacity')
								Cap = A.Capacity;
							else
								Cap = [];
							end
							constructByAdjacent(Adj, Cap);
						else
							error('error:[%s] unknown data structure.', calledby);
						end
					end
				case 2
					if istable(A)
						constructByEdgeNodeTable(A, C);
					else
						constructByAdjacent(A, C);
					end
				otherwise
			end
			
			if ~isempty(this.Head) && isempty(this.link_id)
				this.link_id = spalloc(this.NumberNodes, this.NumberNodes, this.NumberEdges);
				for eid=1:length(this.Head)
					% link e from head(e) to tail(e)
					this.link_id(this.Head(eid),this.Tail(eid))=eid;
				end
			end
			
			function constructByEdgeNodeTable(edge_table, node_table)
				this.Head = edge_table.EndNodes(:,1);
				this.Tail = edge_table.EndNodes(:,2);
				num_links = height(edge_table);
				if contains('Weight', edge_table.Properties.VariableNames)
					this.LinkWeight = edge_table.Weight;
				else
					this.LinkWeight = ones(num_links,1);
				end
				if isempty(node_table)
					num_nodes = max(max(edge_table.EndNodes));
				else
					num_nodes = height(node_table);
				end
				this.Adjacent = spalloc(num_nodes, num_nodes, num_links);
				for e = 1:num_links
					s = this.Head(e);
					t = this.Tail(e);
					this.Adjacent(s,t) = this.LinkWeight(e);
				end
				if contains('Capacity', edge_table.Properties.VariableNames)
					this.Capacity = spalloc(num_nodes, num_nodes, num_links);
					for e = 1:num_links
						s = this.Head(e);
						t = this.Tail(e);
						this.Capacity(s,t) = edge_table.Capacity(e);
					end
				end
			end
			
			function constructByAdjacent(Adj, Cap)
				this.Adjacent = Adj;
				[this.Head, this.Tail, this.LinkWeight] = find(this.Adjacent);
				if ~isempty(Cap)
					this.Capacity = Cap;
				end
			end
		end
		%% Shallow copy only
		
		%% Add new edges to the graph
		% TODO: rename.
		%
		% After adding new edges (new_head,new_tools), the colum-indexing
		% rule of edges might be broken, since the new edges might precede
		% the old edges by the column-indexing rules in the adjacency matrix.
		% For convienence, we keep the indices of old edges not changed
		% (which have been recorded as |link_id| in the object), while
		% appending new edges with new indices.
		%
		%   this = Update(this, head, tail, props)
		function this = Update(this, head, tail, props)
			if isempty(head) || isempty(tail)
				return;
			end
			prev_num_edges = this.NumberEdges;
			this.Head = [this.Head; head];
			this.Tail = [this.Tail; tail];
			[~, b] = unique([this.Head this.Tail], 'rows');
			if length(b) < length(this.Head)
				error('error: duplicate edges.');
			end
			this.LinkWeight = [this.LinkWeight; props.Weight];
			num_nodes = max([this.Head; this.Tail]);
			if size(this.Adjacent) < num_nodes  % expansion matrix
				this.Adjacent(num_nodes, num_nodes) = 0;
				this.link_id(num_nodes, num_nodes) = 0;
			end
			idx = prev_num_edges + 1;
			for i=1:length(head)
				if this.link_id(head(i),tail(i)) == 0
					this.link_id(head(i),tail(i))= idx;
					idx = idx + 1;
				else
					warning('edge (%d,%d) exists, update properties.', head(i), tail(i));
				end
				this.Adjacent(head(i),tail(i)) = props.Weight(i);
			end
			
			if isfield(props, 'Capacity')
				if size(this.Adjacent) < num_nodes  % expansion matrix
					this.Capacity(num_nodes, num_nodes) = 0;
				end
				for i = 1:length(head)
					this.Capacity(head(i),tail(i)) = props.Capacity(i);
				end
			elseif ~isempty(this.Capacity)
				error('error: property <Capacity> should be specified.')
			end
		end
		
		function this = Replace(this, head, tail, props)
			this.Head = head;
			this.Tail = tail;
			this.LinkWeight = props.Weight;
			num_nodes = max([this.Head;this.Tail]);
			num_edges = length(head);
			this.Adjacent = spalloc(num_nodes, num_nodes, num_edges);
			this.link_id = spalloc(num_nodes, num_nodes, num_edges);
			for eid=1:num_edges
				this.link_id(head(eid),tail(eid))=eid;
				this.Adjacent(head(eid),tail(eid)) = props.Weight(eid);
			end
			if isfield(props, 'Capacity')
				this.Capacity = spalloc(num_nodes, num_nodes, num_edges);
				for eid = 1:num_edges
					this.Capacity(head(eid),tail(eid)) = props.Capacity(eid);
				end
			elseif ~isempty(this.Capacity)
				error('error: property <Capacity> should be specified.')
			end
		end
		
		%%
		% Remove links and reallocate link index.
		function [rm_nodes, rm_links] = Remove(this, b_rm_nodes, b_rm_links)
			if ~islogical(b_rm_links)
				temp = false(this.NumberLinks,1);
				temp(b_rm_links) = true;
				b_rm_links = temp;
			end
			if ~isempty(b_rm_nodes)
				if ~islogical(b_rm_nodes)
					temp = false(this.NumberNodes,1);
					temp(b_rm_nodes) = true;
					b_rm_nodes = temp;
				end
				for i = 1:length(b_rm_nodes)
					if b_rm_nodes(i)
						b_rm_links = b_rm_links | (this.Head == i | this.Tail == i);
					end
				end
			end
			link_rm_index = find(b_rm_links);
			if isempty(link_rm_index)
				rm_nodes = [];
				if nargout >= 2
					rm_links = [];
				end
				return;
			end
			for eid=(link_rm_index(:))'
				h = this.Head(eid);
				t = this.Tail(eid);
				this.Adjacent(h,t) = 0;
			end
			b_rm_nodes = sum(this.Adjacent,1)==0 & (sum(this.Adjacent,2)==0)';
			this.Adjacent(b_rm_nodes, b_rm_nodes) = [];
			this.Head(b_rm_links) = [];
			this.Tail(b_rm_links) = [];
			this.link_id = spalloc(this.NumberNodes, this.NumberNodes, this.NumberEdges);
			for eid=1:this.NumberEdges
				this.link_id(this.Head(eid),this.Tail(eid))=eid;
			end
			rm_nodes = b_rm_nodes;
			if nargout >= 2
				rm_links = b_rm_links;
			end
			if ~isempty(this.Capacity)
				this.Capacity(b_rm_nodes, b_rm_nodes) = [];
			end
		end
	end
	
	%% Properties
	methods
		function n = get.NumberNodes(this)
			n = size(this.Adjacent,1);
		end
		
		function m = get.NumberEdges(this)
			m = length(this.Head);
		end
		
		% 		function tbl = get.EdgeTable(this)
		% 			[head, tail, idx] = find(this.link_id');
		% 			tbl = table(idx, [tail,head], 'VariableName', {'Index', 'EndNodes'});
		% 		end
		
		%% Set Weight
		% TODO: rename as SetWeight.
		function SetAdjacent(this, i, j, w)
			% check data vailidity.
			if this.Adjacent(i,j) == 0
				error('error: <%d,%d> does not exist.', i, j);
			end
			this.Adjacent(i,j) = w;
		end
		
		% (Deprecated, now directly access from the property name.)
		%         function w = GetAdjacent(this, i, j)
		%             w = this.Adjacent(i,j);
		%         end
		%% Incident Matrix
		% flow-out is negative, flow-in is positive
		function sp_mat = GetIncidentMatrix(this, bFull)
			L = this.NumberEdges;
			n = [this.Head; this.Tail];
			e = [(1:L)'; (1:L)'];
			b = [-ones(L,1); ones(L,1)];
			sp_mat = sparse(n, e, b);
			if nargin >= 2 && bFull
				sp_mat = full(sp_mat);
			end
		end
	end
	
	methods
		%% Edge Indexing
		% idx = IndexEdge(this, src, dest)
		% get the link index by the specified link head and tail. |src| and |dest| are
		% scalar or vector, represent the link heads and tails.
		%
		% Example:
		%		[src dest] = IndexEdge(this)
		%		[src dest] = IndexEdge(this, id)
		%							   Get the link head and tail by the speified link |id|,
		%							   |id| is a vector.
		%		idx = IndexEdge(this, s, t)
		function [arg1_out,arg2_out] = IndexEdge(this, arg1_in, arg2_in)
			if nargin == 1
				arg1_out = this.Head;
				arg2_out = this.Tail;
			elseif nargin == 2
				idx = arg1_in;
				arg1_out = this.Head(idx);
				arg2_out = this.Tail(idx);
			elseif nargin==3
				s = arg1_in;
				t = arg2_in;
				if isempty(s) || isempty(t)
					arg1_out = double.empty(0,1);
					return;
				end
				if isscalar(s) && ~isscalar(t)
					s = s*ones(size(t));
				end
				if isscalar(t) && ~isscalar(s)
					t = t*ones(size(s));
				end
				if length(s)~=length(t)
					error('arguments src and dest are vectors, which must have the same length');
				end
				arg1_out = zeros(length(s),1);
				for i=1:length(s)
					arg1_out(i) = this.link_id(s(i),t(i));
				end
			else
				error('invalid number of input arguments');
			end
		end
		
		%% Select Candidate Paths
		% This function finds the shortest paths from single source node to single
		% destination node. It also allows to find the shortest path from a single source
		% to a set of nodes, which is useful when multi-source service is considered.
		%
		%   [path_list, k] = CandidatePaths(this, K, src, dest_set)
		%
		% *Input Arguments*
		%
		% |K|: Number of candidate paths to be selected.
		%
		% |src|: source node of the candidate paths.
		%
		% |options|: options for select candidate paths, including:
		%		(1) _delay_: delay constraint for path,
		%		(2) _delay_opt_: if this option is set as |'BandwidthDependent'|,
		%				link delay depends on residual bandiwdth,
		%		(3) _mid_set_: the middle nodes, one of which must be traversed by
		%				the candidate path.
		%
		% *Output Arguments*
		%
		% |path_list|: List of candidate paths. If the source and destination are not
		% connected, or no path meet the delay requirement, |path_list| is returned with
		% an empty array.
		%
		% |k|: actual number of the selected paths, it may be less than K.
		function [path_list, k] = CandidatePaths(this, K, src, dest_set, options)
			global DEBUG;
			defaultopts = struct(...
				'DelayConstraint', inf,...
				'DelayModel', LinkDelayOption.BandwidthPropotion,...
				'MiddleNodes', []);
			if nargin <= 5
				options = defaultopts;
			else
				options = structupdate(defaultopts, options);
			end
			if options.DelayConstraint == 0
				options.DelayConstraint = inf;
			end
			
			k = 0;
			graph = this.copy;  % for all flows, the graph should be the same.
			path_list = cell(K,1);
			while k < K
				%%%
				% If the shortest path does not traverse the middle nodes, then we need to
				% find a shortest path from |src| to |mid_set|, then to |dest|. This
				% search procedure may take more time. For example, if there are |n|
				% middle nodes in |mid_set|, and |m| nodes in |dest_set|, there are
				% $m\times n$ combinations of the shortest path.
				% We perform the search in the following way. First we run shortest path
				% algorithm to find if the shortest path traverse one of the middle nodes.
				% If true, the procedure is fininshed. Otherwise, we need to cintinue the
				% search procedure, by enumerating all the shortest path from |src| to
				% middle nodes and middle nodes to |dest_set|. thus, the shortest path
				% that traverse at least a mniddle nodes can be found.
				%
				% We use dijkstra's algorithm to calculate the shortest paths instead of
				% All Pair Shortest Path algorithm (_i.e._ Floyd's algorithm). Considering
				% |n| middle nodes and |m| dest nodes, then dijkstra's algorithm need to
				% be executed for |n+1| times, so the complexity is
				% $\mathcal{O}(N^2(n+1))$(un-optimized), where N is the total number of
				% nodes in the network. On the other hand, the complexity of Floyd's
				% algorithm is $\Theta(N^3)$. Therefore, we choose  Dijkstra's algorithm.
				% Besides, when mid_set is small, the complexity is much lower; and when
				% mide_set is large, there is huge chance that we first find a shortest
				% path from source to dest_set that travese the mid_set, and there is no
				% need to perform the enumeration.
				
				% TODO, the reverse direction
				% Note: ShortestPath may find the graph is not all connected if k > 1,
				% since some links may have been removed from the graph.
				[pn, flag, path_delay] = graph.SingleSourceShortestPaths(src,dest_set);
				if flag == 0 && ~isempty(options.MiddleNodes)
					if isempty(intersect(options.MiddleNodes, pn))
						[pn, flag, path_delay] = ...
							graph.ForcedShortestPath(src, dest_set, options.MiddleNodes);
					end
				end
				if flag ~= 0
					if k==0
						if ~isempty(DEBUG) && DEBUG
							warning('graph does not connected between %d and %d.', ...
								src, dest_set(1));
						end
						path_list = [];
						return;
					else
						break;
					end
				else
					if path_delay > options.DelayConstraint
						if k == 0
							warning('no path meet the delay constraint.');
							path_list = [];
							return;
						else
							break;
						end
					end
					bw = inf;
					for j=1:(length(pn)-1)
						bw = min(bw, graph.Capacity(pn(j),pn(j+1)));
					end
					for j=1:(length(pn)-1)
						graph.Capacity(pn(j),pn(j+1)) = graph.Capacity(pn(j),pn(j+1)) - bw;
						if options.DelayModel == LinkDelayOption.BandwidthInverse
							% link weight is inverse propotional to the available bandwidth
							graph.Adjacent(pn(j),pn(j+1)) = ...
								PhysicalNetwork.LinkDelay(options.DelayModel, graph.Capacity(pn(j),pn(j+1)));
						else
							% link weight is not changed untill the capacity is depleted.
							if graph.Capacity(pn(j),pn(j+1)) <= 0
								graph.Adjacent(pn(j),pn(j+1)) = 0;
							end
						end
					end
					k = k + 1;
					path_list{k} = Path(pn, bw, path_delay);
					% Note: the initialization of bandwidth has no effect.
				end
			end
			path_list(k+1:end) = [];
			path_list = PathList(path_list);
		end
		
		%%%
		% * *ForcedShortestPath*
		% the shortest path must pass one of the node in the |mid_set|.
		function [path, flag, cost] = ...
				ForcedShortestPath(this, source, targets, mid_set, options)
			path = [];
			cost = [];
			options.b_exhausted = true;
			[path_list1, flag, path_delay] = ...
				this.SingleSourceShortestPaths(source, mid_set, options);
			if flag ~= 0
				return;
			end
			
			% filter middle nodes
			l1 = path_delay<inf;
			for i = find(l1)
				if find(mid_set == path_list1{i}(1:(end-1)),1)
					l1(i) = false;
				end
			end
			mid_set = mid_set(l1);
			path_list1 = path_list1(l1);
			path_delay = path_delay(l1);
			%
			cost_table = zeros(nnz(l1), length(targets));
			path_table = cell(nnz(l1), length(targets));
			for i = 1:length(mid_set)
				[path_table(i,:), ~, cost_table(i,:)] = ...
					this.SingleSourceShortestPaths(mid_set(i), targets, options);
			end
			dist = min(cost_table);
			if dist == inf
				flag = -1;
				return;
			end
			[mid_src_index, dest_index] = find(cost_table==dist,1);
			flag = 0;
			path = [path_list1{mid_src_index}, path_table{mid_src_index,dest_index}(2:end)];
			cost = path_delay(mid_src_index) + cost_table(mid_src_index,dest_index);
		end
		
		%% Single Source Shortest Paths
		% |options|: |direction|,  |b_exhausted, algorithm|
		%
		% |path|: |path| is a row vector. when only one path is computed, |path| is a
		%       vector of nodes on the shortest path. If a shortest path tree is computed,
		%       |path| is a cell vector contains paths from |source| to each target node.
		% |flag|: if find a shortest path or build a shortest path tree (partially),
		% |flag=0|. Otherwise, |flag=-1|.
		function [path, flag, cost] = ...
				SingleSourceShortestPaths(this, source, targets, options)
			if nargin <= 3 || ~isfield(options, 'algorithm')
				options.algorithm = 'dijkstra';
			end
			if nargin <= 3 || ~isfield(options, 'direction')
				options.direction = 0;
			end
			if nargin <= 3 || ~isfield(options, 'b_exhausted')
				options.b_exhausted= false;
			end
			if nargin <=2 || isempty(targets)
				targets = 1:this.NumberNodes;
				targets(source)=[];
				options.b_exhausted = true;
			elseif ~isrow(targets)
				targets = reshape(targets, 1, numel(targets));
			end
			idx = find(targets==source);
			if ~isempty(idx)
				warning('destination must be different from source.');
				targets(idx)=[];
			end
			if max([source targets])>this.NumberNodes || min([source targets]) < 1
				error('argument error: invalid node index.');
			end
			
			switch options.algorithm
				case 'dijkstra'
					[cost, prev, flag] = this.algorithmDijkstra(source, targets, options);
				case 'bellman-ford'
					[cost, prev, flag] = this.algorithmBellmanFord(source, targets);
				otherwise
					error('error: invalid algorithm specification (%s).', ...
						options.algorithm);
			end
			
			%% TODO (back) trace the shortest path
			if ~options.b_exhausted
				if flag < 0
					path = [];
					flag = -1;
					return;
				end
				targets = flag;
			else
				if isempty(find(cost<inf,1))
					path = {[]};
					flag = -1;
					return;
				end
			end
			flag = 0;
			path_list = cell(length(targets),1);
			for i = 1:length(targets)   % i~=source
				if cost(i) == inf
					continue;
				end
				p = zeros(1,this.NumberNodes);
				p(1) = targets(i);
				j = 1;
				while p(j) ~= source
					node = prev(p(j));
					j = j+1;
					p(j) = node;
				end
				if options.direction
					% Since the direction is dest to src, backtrace is not required.
					path_list(i) = {p(1:j)};
				else
					path_list(i) = {p(j:-1:1)};
				end
			end
			if ~options.b_exhausted
				path = path_list{1};
			else
				path = path_list;
			end
		end
		
		%%%
		% Floyd's algorithm
		% |path|: |path(i,j)| records the next hop of the shortest path between |i| and
		%         |j|.
		% |cost|: records the cost of the shortest path between |i| and |j|.
		% |flag|: when |flag=false|, the graph is not all-connected. If |i| and |j| is not
		% connected, |path(i,j)=-1|
		function [cost, flag, path] = AllPairShortestPath(this, options)
			N = this.NumberNodes;
			this.mat_path = ones(N)*(-1);
			this.mat_cost = this.Adjacent;
			for i = 1:N
				for j = 1:N
					if this.mat_cost(i,j) == 0
						this.mat_path(i,j) = 0;
					elseif this.mat_cost(i,j) < inf
						this.mat_path(i,j) = i;      % initialize next hop
					end
				end
			end
			
			for k = 1:N
				for i = 1:N
					for j = 1:N
						new_cost = this.mat_cost(i,k)+this.mat_cost(k,j);
						if this.mat_cost(i,j) > new_cost
							this.mat_cost(i,j) = new_cost;
							this.mat_path(i,j) = this.mat_path(k,j);
						end
					end
				end
				if nargin>=2 && isfield(options, 'Display') && strcmpi(options.Display, 'on')
					disp('Info: cost matrix:');
					disp(this.mat_cost);
					disp('Info: last hop matrix:');
					disp(this.mat_path);
				end
			end
			if nargin >=2
				if find(this.mat_cost==inf)
					flag = true;
				else
					flag = false;
				end
			end
			
			this.mat_cpath = cell(N);
			for i = 1:N
				for j = 1:N
					this.mat_cpath{i,j} = zeros(N-1,1);
					lp = 1;
					this.mat_cpath{i,j}(1) = j;
					last_hop = this.mat_path(i,j);
					if last_hop == -1
						% (i,j) is not connected.
						this.mat_cpath{i,j} = [];
					end
					while last_hop<=0  % have not reach the source
						lp = lp + 1;
						this.mat_cpath{i,j}(lp) = last_hop;
						last_hop = this.mat_path(i,last_hop);
					end
					this.mat_cpath{i,j} = this.mat_cpath{i,j}(lp:-1:1);
					if nargin>=2 && isfield(options, 'Display') && strcmpi(options.Display, 'on')
						fprintf('Info: shortest path from %d to %d:\n', i,j);
						disp(this.mat_cpath{i,j}');
					end
				end
			end
			
			if nargin >= 1
				cost = this.mat_cost;
			end
			if nargin>=3
				if isfield(options, 'ExplicitPath') && options.ExplicitPath
					path = this.mat_cpath;
				else
					path = this.mat_path;
				end
			end
			% *FIXME* only works for non-negative weighted graph.
			%             path2 = ones(N)*(-1);
			%             for i = 1:N
			%                 for idj = 1:N
			%                     if i==j
			%                         path2(i,j) = 0;
			%                         continue;
			%                     end
			%                     [~, neighbor_j] = find(this.Adjacent(:,j)'<inf);
			%                     %%%
			%                     % Since
			%                     % $d_{i,j}^{(k)}=\min{d_{i,j}^{(k-1)},d_{i,k}^{(k-1)}+d_{k,j}^{(k-1)}}$
			%                     % the maximum |k| that satisfies the equation, must be the last hop of
			%                     % pair $(i,j)$.
			%                     for k = neighbor_j
			%                         if k~=j && cost(i,j) == cost(i,k) + this.Adjacent(k,j)
			%                             path2(i,j) = k;
			%                             break;
			%                         end
			%                     end
			%                 end
			%             end
		end
		
		function plot(this, line_spec)
			W = this.Adjacent;
			W(W==inf) = 0;
			g = digraph(W);
			[~,~,edge_lables]= find((this.link_id)');
			if nargin >= 2
				g.plot(line_spec, 'EdgeLabel', edge_lables);
			else
				g.plot('EdgeLabel', edge_lables);
			end
		end
	end
	
	methods (Access=private)
		%% Find Shortest Path by Dijkstra's Algorithm
		% _algorithmDijkstra_ finds a shortest path (tree) from single _source_ to
		% _targets_, by using the graph's adjacent/weight matrix. If there is no direct
		% connection between two nodes in the grpah, the value in the adjacent matrix is
		% |inf|(|0| can be automatically converted to |inf| when calculate the shortest
		% path by _DirectedGraph.d_).
		%
		%   [path, flag, cost] = algorithmDijkstra(source, targets, options)
		%
		% *Input Arguments*
		%
		% |source|: the single source node, which is the root of the shortest path tree;
		%
		% |targets|: target nodes;
		%
		% # If |targets| is not specified, this method build the shortest path tree from
		%     the source to all other nodes.
		% # If |targets| is a vector and the option |b_exhausted = false|, this function
		% return the path from source to the closet nodes in |targets|.
		%
		% |options|: a structure including the following fields:
		%
		% # |direction|: if this field is specified and |direction=1|, the algorithm find
		%       the shortest path from |targets| to |source|. Otherwise, by default, it
		%       finds the shortest path form |source| to |targets|.
		% # |b_exhausted|: if this field is specified and |b_exhausted=true|, the
		%       algorithm will build the shortest path tree from |source| to |targets|.
		%       Otherwise, it only find the path from source to the closet nodes in
		%       |targets|.
		%
		% * *Output Arguments*
		%
		% |cost|: when only one path is computed, |cost| is a scalar for the cost of the
		%       shortest path. If a shortest path tree is calculated, |cost| is a vector
		%       for the cost from |source| to each target node.
		%
		% |prev|: last hop on the shortest path of the visited nodes.
		%
		% |flag|: if a shortest path is found, |flag=target|; if a shortest path tree is
		%       found, |flag=0|; Otherwise, _source_ and _targets_ are not connected,
		%       |flag=-target|, where |target| records the first unconnected node.
		%
		% *NOTE*: (1) if the graph is not all-connected, this method return an error flag
		%       instead of throwing an error. (2) this is a private function, all
		%       arguments must be provided.
		function [cost, prev, flag] = algorithmDijkstra(this, source, targets, options)
			% assertion
			N = this.NumberNodes;		% the number of nodes in the graph
			if options.direction==0
				W = DirectedGraph.d(this.Adjacent);
			else
				% reverse the direction of edges in the graph
				W = DirectedGraph.d(this.Adjacent');
			end
			
			dist = full(W(source,:));   % distance from source node to other joint nodes
			b_permanent = zeros(1,N);	% permanent label for nodes
			dist(source) = 0;
			b_permanent(source) = 1;
			prev = ones(N,1)*source;    % the recode of previous node on the shortest path
			if options.b_exhausted
				count = length(targets);
				b_dest = zeros(1,N);
				b_dest(targets) = 1;
			end
			flag = 0;
			
			for i = 1:N-1
				temp_dist = dist;
				temp_dist(b_permanent==1) = inf;
				[min_dist,index] = min(temp_dist);
				if min_dist == inf
					% the temporary nodes include the targets are not reachable from the
					% source.
					flag = -index;
					break;
				end
				b_permanent(index) = 1;
				if ~isempty(find(index == targets,1)) && ~options.b_exhausted
					% when the destination is reached, the search procedure is terminated.
					flag = index;
					break;
				end
				for k = find(W(index,:)<Inf)
					new_dist = dist(index)+W(index,k);
					if dist(k) > new_dist
						dist(k) = new_dist;
						prev(k) = index;
					end
				end
				
				if options.b_exhausted && b_dest(index)==1
					count = count - 1;
					if count == 0
						break;
					end
				end
			end
			
			if options.b_exhausted
				cost = dist(targets);   % return the distance to the target nodes
			else
				cost = dist(index);     % return the distance to the closet target node
			end
		end
		
		function [dist, path, flag] = algorithmBellmanFord(this, source)
			N = this.NumberNodes;
			L = this.NumberEdges;
			W = DirectedGraph.d(this.Adjacent);
			dist = ones(N,1)*inf;
			dist(source) = 0;
			path = zeros(N,1);
			
			for i=1:(N-1)
				for j=1:L
					u = this.Head(j);
					v = this.Tail(j);
					if dist(u) > dist(v) + W(v,u)
						dist(u) = dist(v) + W(v,u);
						path(u) = v;    % predecessor
					end
				end
			end
			% determine if negative cycle exists.
			flag = true;
			for j=1:L
				u = this.Head(j);
				v = this.Tail(j);
				if dist(u) > dist(v) + W(v,u)
					flag = false;
					break;
				end
			end
		end
	end
	
	methods(Access=protected, Static)
		function l = d(dist)
			dist(dist==0)=inf;
			l = dist;
		end
	end
end

