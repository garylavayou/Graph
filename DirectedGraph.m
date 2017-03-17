%% Directed Graph
% directed graph

%%
classdef DirectedGraph < matlab.mixin.Copyable
    properties
        Adjacent;       % matrix link weight
        Capacity;       % matrix link capacity
        
        Head;
        Tail;
        LinkWeight;
%         LinkCapacity;
        
        NumberNodes;
        NumberEdges;
        
    end

    properties (Access = private)
        link_id;         % inverse index of links
    end
    
    %% Constructor
    methods
        function this = DirectedGraph(A, C)
            switch nargin
                case 1
                    if isa(A, 'DirectedGraph')
                        this = A.copy;
                    elseif isa(A, 'digraph')
                        this.Adjacent = A.adjacency();
                        %% Edges Index Inconsistent
                        % In the digraph, edges are indexed by rows, while matlab matrix
                        % is indexed by column. Since in Directed Graph, we use matlab
                        % matrix to save the adjacent, the edge index of digraph cannot be
                        % directly used, i.e., 
                        %   this.Head = A.findnode(A.Edges.EndNodes(:,1));
                        %   this.Tail = A.findnode(A.Edges.EndNodes(:,2));
                        % Instead, we use adjacent matrix column index as follows.
                        [this.Head, this.Tail] = find(this.Adjacent~=0);
                        this.LinkWeight = zeros(size(this.Head));
                        this.Capacity = spalloc(A.numnodes,A.numnodes,A.numedges);
                        for e = 1:length(this.Head)
                            s = this.Head(e);
                            t = this.Tail(e);
                            eid = A.findedge(s,t);
                            this.Capacity(s,t) = A.Edges.Capacity(eid);
                            this.Adjacent(s,t) = A.Edges.Weight(eid);
                            this.LinkWeight(e) = A.Edges.Weight(eid);
                        end
                    elseif isnumeric(A)
                        this.Adjacent = A;
                        [this.Head, this.Tail, this.LinkWeight] = find(this.Adjacent);
                        this.Capacity = double(this.Adjacent~=0);
                    end
                case 2
                    this.Adjacent = A;
                    [this.Head, this.Tail, this.LinkWeight] = find(this.Adjacent);
                    this.Capacity = C;
                otherwise
            end
            
            if ~isempty(this.Head) && isempty(this.link_id)
                this.link_id = spalloc(this.NumberNodes, this.NumberNodes, this.NumberEdges);
                for e=1:length(this.Head)
                    % link e from head(e) to tail(e)
                    this.link_id(this.Head(e),this.Tail(e))=e;
                end
            end
        end
        %% Shallow copy
    end
    
    %% Properties
    methods
        function n = get.NumberNodes(this)
            n = size(this.Adjacent,1);
        end
        
        function m = get.NumberEdges(this)
            m = length(this.Head);
        end
        
        function SetAdjacent(this, i, j, w)
            this.Adjacent(i,j) = w;
            %% TODO
        end
        
        function w = GetAdjacent(this, i, j)
            w = this.Adjacent(i,j);
        end
    end
    
    methods
        %% Edge Indexing
        % idx = IndexEdge(this, src, dest)
        % get the link index by the specified link head and tail. |src| and |dest| are
        % scalar or vector, represent the link heads and tails.
        % 
        % [src dest] = IndexEdge(this, id)
        % [src dest] = IndexEdge(this)
        % get the link head and tail by the speified link |id|, |id| is a vector.
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
                if isempty(s)
                    s = 1:this.NumberNodes;
                end
                if isempty(t)
                    t = 1:this.NumberNodes;
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
        % |options|: options for select candidate paths. includinf _delay_ delay
        % constraint for path, and _delay_opt_ 'BandwidthDependent': delay dependence on
        % bandiwdth. 
        %
        % *Output Arguments*
        %
        % |path_list|: List of candidate paths. If the source and destination are not
        % connected, or no path meet the delay requirement, |path_list| is returned with
        % an empty array.
        %
        % |k|: actual number of the selected paths, it may be less than K.
        function [path_list, k] = CandidatePaths(this, K, src, dest_set, options)
            if nargin < 5 
                options.delay = inf;
                options.delay_opt = LinkDelayOption.BandwidthPropotion;
            else
                if isempty(options.delay) || options.delay == 0
                    options.delay = inf;
                end
                if isempty(options.delay_opt)
                    options.delay_opt = LinkDelayOption.BandwidthPropotion;
                end
            end
            
            k = 0;
            graph = DirectedGraph(this);
            path_list = cell(K,1);
            while k < K
                % TODO, the reverse direction
                % Note: ShortestPath may find the graph is not all connected if k > 1, since
                % some links may have been removed from the graph.
                [pn, flag, path_delay] = graph.ShortestPath(src,dest_set); 
                if flag ~= 0
                    if k==0
                        warning('graph does not connected between %d and %d.', src, flag);
                        path_list = [];
                        return;
                    else
                        break;
                    end
                else
                    if path_delay > options.delay
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
                        if options.delay_opt == LinkDelayOption.BandwidthInverse
                            % link weight is inverse propotional to the available bandwidth
                            graph.Adjacent(pn(j),pn(j+1)) = ...
                                PhysicalNetwork.LinkDelay(options.delay_opt, graph.Capacity(pn(j),pn(j+1)));
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
        
        %% Shortest Path
        % _ShortestPath_ finds a shortest path from "source" to "destination". The matrics
        % is stored as the graph's adjacent/weight matrix. If two nodes are not connected,
        % the value in the adjacent matrix is 0/inf (0 can be automatically converted to
        % inf when calculate the shortest path). 
        %
        %   [path, flag] = ShortestPath(source, targets, dir)
        %
        % *Input Arguments*
        %
        % |source|: Source node, i.e. the root of the shortest path tree;
        %
        % |targets|: destination nodes;
        %
        % # If |targets| is not specified, this method find the shortest path to
        %     most distant node.
        % # If |targets| is a vector, this function return the path from source to
        %     the closet nodes in |targets|.
        %
        % |dir|: whether the path is from |source| to |targets| (|dir=0|) or from
        %        |targets| to |source|(|dir=1|).
        %
        % * Output Arguments*
        %
        % |path|: a path consists of a sequence of nodes from source to destination.
        %
        % |flag|: 0 - find a shortest path;
        %       otherwise, the source and destination are not connected, flag record
        %       the first unconnected node.
        %
        % |cost|: cost of shortest path.
        %
        % *NOTE1*: Difference from _ShortestPathTree_: (1) this method retruns a path not
        % a set of path, (2) if the graph is not all-connected, this method return a error
        % flag instead of throwing an error.
        function [path, flag, cost] = ShortestPath(this, source, targets, dir)
            % assertion
            N = this.NumberNodes;		% the number of nodes in the graph
            if nargin <=2 || isempty(targets)
                targets = 1:N;
                targets(source)=[];
            elseif iscolumn(targets)
                targets = reshape(targets, 1, numel(targets));
            end
            if nargin <= 3
                dir = 0;
            end
            
            idx = find(targets==source);
            if ~isempty(idx)
                warning('destination must be different from source.');
                targets(idx)=[];
            end
            if max(targets)>N || source > N || min(targets) < 1 || source < 1
                error('argument error: invalid node index.');
            end
            if dir==0
                W = this.Adjacent;
            else
                W = this.Adjacent';     % reverse the direction of edges in the graph
            end
            
            
            dist = DirectedGraph.d(full(W(source,:))); % distance from source node to other joint nodes
            b_permanent = zeros(1,N);	% permanent label for nodes
            dist(source) = 0;
            b_permanent(source) = 1;
            prev = ones(N,1)*source;    % the recode of previous node on the shortest path
            
            for i = 1:N-1
                temp_dist = dist;
                temp_dist(b_permanent==1) = inf;
                [min_dist,index] = min(temp_dist);
                if min_dist == inf
                    % the temporary nodes include the destination are not reachable to the
                    % source.
                    flag = index;
                    path = [];
                    cost = [];
                    return;
                end
                b_permanent(index) = 1;
                if ~isempty(find(index == targets,1))
                    % when the destination is reached, the search procedure is terminated.
                    flag = 0;
                    dest = index;
                    break;
                end
                for k = find(DirectedGraph.d(W(index,:))<Inf)
                    new_dist = dist(index)+DirectedGraph.d(W(index,k));
                    if dist(k) > new_dist
                        dist(k) = new_dist;
                        prev(k) = index;
                    end
                end
            end
            
            % (back) trace the shortest path
            path = zeros(1,N);
            path(1) = dest;
            j = 1;
            while path(j) ~= source
                node = prev(path(j));
                j = j+1;
                path(j) = node;
            end
            if dir == 0
                path = path(j:-1:1);
            else
                path = path(1:j);
            end
            cost = dist(targets);
        end
        
        function plot(this)
            g = digraph(this.Adjacent);
            g.plot;
        end
    end
        
    methods(Access=private, Static)
        function l = d(dist)
            dist(dist==0)=inf;
            l = dist;
        end
    end
end

