%% Graph

%%
%   Graph < handle
classdef Graph < handle
    %GRAPH Summary of this class goes here
    %   Detailed explanation goes here
    properties(Constant=true)
        THETA = 0.5;
    end
    properties
        Adjacent;   % capacity of edges: if edge(i,j) exists, Adjacent(i,j) is initialized to 1, otherwise 0.

        Head;
        Tail;
        Incidence;
        Inverse;    % inverse edge index (i,j) -> e
        
        Capacity;   % capacity of edges;
        DistAdjacent;
        
        Size;       % Read only
        EdgeNumber;
    end
    
    methods
        function this = Graph(graph)
            switch nargin
                case 1
                    if isa(A, 'Graph')
                        this.Copy(graph);
                    elseif isnumeric(A)
                        this.Adjacent = A;
                        this.Renew;
                    end
                otherwise
            end
        end
        
        function n = get.Size(this)
            n = size(this.Adjacent,1);
        end
        
        function n = get.EdgeNumber(this)
            n = length(this.Capacity);
        end
        
        function IncreaseCapacity(this, a)
            this.Adjacent = this.Adjacent*a;
            this.Capacity = this.Capacity*a;
            this.DistAdjacent = this.DistAdjacent/a;
        end
        
        function Copy(this, graph)
            this.Adjacent = graph.Adjacent;
            this.Head = graph.Head;
            this.Tail = graph.Tail;
            this.Incidence = graph.Incidence;
            this.Capacity = graph.Capacity;
            this.Inverse = graph.Inverse;
            this.DistAdjacent = graph.DistAdjacent;
        end
        
        function SetAdjacent(this, Adj)
            this.Adjacent = Adj;
            this.Renew;
        end
        function Renew(this)
            [this.Head, this.Tail, this.Capacity] = find(this.Adjacent);
            this.Inverse = zeros(this.Size);
            for e=1:length(this.Capacity)
                % link e from head(e) to tail(e)
                this.Inverse(this.Head(e),this.Tail(e))=e;
            end
            this.DistAdjacent = 1./this.Adjacent;
            for i = 1:this.Size
                this.DistAdjacent(i,i) = 0;
            end
        end
        %% Select Candidate Paths
        % This function allows select paths from single source node to single destination
        % node. It also allows select paths from single source to multiple destinations if
        % mutiple target nodes are specified.
        %
        %   [path_list, k] = CandidatePaths(this, K, src, dest_set)
        %
        % *Input Arguments*
        %
        % |K|: Number of candidate paths to be selected.
        %
        % |src|: source node of the candidate paths.
        %
        % *Output Arguments*
        %
        % |path_list|: List of candidate paths.
        %
        % |k|: actual number of the selected paths, it may be less than K.
        function [path_list, k] = CandidatePaths(this, K, src, dest_set)
            k = 0;
            graph = Graph;
            graph.Copy(this);
            path_list = cell(K,1);
            while k < K
                [pn, flag] = graph.ShortestPath(src,dest_set); % TODO, the reverse direction
                if flag ~= 0
                    if k==0
                        error('error: graph does not connected between %d and %d.', src, flag);
                    else
                        break;
                    end
                else
                    bw = inf;
                    for j=1:(length(pn)-1)
                        bw = min(bw, graph.Adjacent(pn(j),pn(j+1)));
                    end
                    for j=1:(length(pn)-1)
                        graph.Adjacent(pn(j),pn(j+1)) = graph.Adjacent(pn(j),pn(j+1)) - bw;
                        graph.DistAdjacent(pn(j),pn(j+1)) = 1/graph.Adjacent(pn(j),pn(j+1));
                    end
                    k = k + 1;
                    path_list{k} = Path(pn, bw);  
                    % Note: the initialization of bandwidth has no effect.
                end
                
            end
        end
        
        %% Shortest Path Tree
        % _ShortestPathTree_ implements _Dijkstra's Algorithm_ to generate a tree that the
        % distance from each node to the root, i.e. the source node (or the reverse
        % direction) is the shortest.
        %
        %   [path_list, distances] = ShortestPathTree(source, targets, dir)
        %
        % *Input Arguments*
        %
        % |source|: Source node, i.e. the root of the shortest path tree;
        %
        % |targets|: A vector represent the set of destination nodes. When all nodes in
        % |targets| has been reached, the spanning of tree is terminated.
        %
        % * If |targets| is empty or not specified. The shortest path tree reaches all
        % nodes in the graph.
        %
        % |dir|: Determine whether the distance is counted from the root to the leaf nodes
        % on the shortest path tree or the reverse direction.
        %
        % * 0 - shortest paths from source to all other nodes(default);
        % * 1 - shortest paths from all other nodes to source;
        %
        % *Output Arguments*
        %
        % |path_list|: The tree is represented as a series of paths between the root and
        % the leaf nodes, organised as a cell vector. Corresponding to the leaf node in
        % |targets|, each cell element is a path consists of a sequence of nodes from
        % root to leaf node.
        %
        % |distances|: the distance of each node in |targets| to the root.
        %
        % *NOTE*: there are more efficient algorithms using Fibonacci heap.
        function [path_list, distances] = ShortestPathTree(this, source, targets, dir)
            % assertion
            N = this.Size;
            if nargin <=2 || isempty(targets)
                targets = 1:N;
                targets(source)=[];
            end
            if nargin <= 3
                dir = 0;
            end
            
            if ~isempty(find(targets==source,1))
                error('argument error: destination must be different from source.');
            end
            if max(targets)>N || source > N || min(targets) < 1 || source < 1
                error('argument error: invalid node index.');
            end
            if dir==0
                W = this.DistAdjacent;
            else
                W = this.DistAdjacent'; % reverse the direction of edges in the graph
            end
            targets = reshape(targets, 1, numel(targets));
            
            %% Dijkstra Algorithm
            % *Initialize*
            %
            % * distances from source node to other nodes. The distance from source |s| to
            % the adjacent node |j| is |W(s,j)|, while the distances from |s| to other
            % nodes is $\infty$.
            dist = W(source,:);        % equal to W(:,source) if dir=1
            %%%
            % * permanent label for nodes. The source node is first labeled as permanent.
            b_permanent = zeros(1,N);
            b_permanent(source) = 1;
            %%%
            % * the record of previous node on the shortest path. This is update during
            % the algorithm, and need not to be initilized.
            prev = ones(N,1)*source;
            count = length(targets);
            b_dest = zeros(1,N);
            b_dest(targets) = 1;
            
            % O(N) = 5N^2 +14N
            for i = 1:N-1
                temp_dist = dist;
                temp_dist(b_permanent==1) = inf;
                %%%
                % *Update*
                % Find a nearest node to source.
                % If a node is permanently labelled, do not consider it.
                [min_dist,index] = min(temp_dist);
                % Here, a connected graph is assumed, otherwise an assertion is needed to
                % break the loop.
                if min_dist == inf
                    error('Error: the graph is not all connected.');
                end
                b_permanent(index) = 1;
                %%%
                % *update the distance function*
                for k = find(W(index,:)<Inf)
                    %%%
                    % Permanent nodes' distance and previous node will not change any more,
                    % so they are excluded from the distance update procedure.
                    %
                    %   if b_permanent(k) == 1
                    %      continue;
                    % end
                    new_dist = dist(index)+W(index,k);         % W(u,u) = 0;
                    %%%
                    % *Note*: if |W(u,u)| is not set to 0, instead it is set to |Inf|,
                    % then |dist(source)=Inf|, thus |dist(source)>t|, and then
                    % |prev(source)=index|. Therefore the result is not correct, which
                    % influence the backtrace procedure if the trace stop condition is
                    % not well set.
                    if dist(k) > new_dist
                        dist(k) = new_dist;
                        prev(k) = index;
                        %%%
                        % For the node added to permanent set just above, its previous node
                        % has been determined in the previous iteration or in the
                        % initialization phase(source node is its previous node).
                    end
                end
                % If the |targets| is always far less than the node set, this assertion is
                % necessary.
                if b_dest(index)==1
                    count = count - 1;
                    if count == 0
                        break;
                    end
                end
            end
            
            %%%
            % *(back) trace the shortest path*
            path_list = cell(length(targets),1);
            k = 1;
            for i = targets   % i~=source
                path = zeros(1,N);
                path(1) = i;
                j = 1;
                while path(j) ~= source
                    node = prev(path(j));
                    j = j+1;
                    path(j) = node;
                end
                if dir
                    %  Since the direction is dest to src, backtrace is not required.
                    path_list(k) = {path(1:j)};
                else
                    path_list(k) = {path(j:-1:1)};
                end
                k = k + 1;
            end
            distances = dist(targets);
        end
        
        %% Shortest Path
        % _ShortestPath_ finds a shortest path from "source" to "destination" (or the
        % reverse direction), by Dijkstra's algorithm.
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
        % |path|: a path consists of a sequence of nodes from source to destination.
        %
        % |flag|: 0 - find a shortest path;
        %       otherwise, the source and destination are not connected, flag record
        %       the first unconnected node.
        %
        %
        % *NOTE1*: Difference from _ShortestPathTree_: (1) this method retruns a path not
        % a set of path, (2) if the graph is not all-connected, this method return a error
        % flag instead of throwing an error.
        function [path, flag] = ShortestPath(this, source, targets, dir)
            % assertion
            N = this.Size;		% the number of nodes in the graph
            if nargin <=2 || isempty(targets)
                targets = 1:N;
                targets(source)=[];
            end
            if nargin <= 3
                dir = 0;
            end
            
            if ~isempty(find(targets==source,1))
                error('argument error: destination must be different from source.');
            end
            if max(targets)>N || source > N || min(targets) < 1 || source < 1
                error('argument error: invalid node index.');
            end
            if dir==0
                W = this.DistAdjacent;
            else
                W = this.DistAdjacent'; % reverse the direction of edges in the graph
            end
            targets = reshape(targets, 1, numel(targets));
            
            dist = W(source,:);     	% distance from source node to other joint nodes
            b_permanent = zeros(1,N);	% permanent label for nodes
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
                    return;
                end
                b_permanent(index) = 1;
                if ~isempty(find(index == targets,1))
                    % when the destination is reached, the search procedure is terminated.
                    flag = 0;
                    dest = index;
                    break;
                end
                for k = find(W(index,:)<Inf)
                    new_dist = dist(index)+W(index,k);
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
            path = path(j:-1:1);
        end
        
        %% Sort node by distance
        % _DistanceOrderedNodes_ returns the sequence of nodes that are sorted by the
        % distance to source node.
        %
        %   [ordered_nodes, distances] = DistanceOrderedNodes(W, source, targets, dir)
        %
        % *Input Arguments* See _Shortest Path Tree_.
        %
        % *Output Arguments*
        %
        % |ordered_nodes|: the nodes in |targets| are sorted in ascending order by
        % distance to the |source|.
        %
        % |distances|: the distances of |ordered_nodes| to the |source|.
        %
        function [ordered_nodes, distances, path_list] = ...
                DistanceOrderedNodes(this, source, targets, dir)
            [path_list, distances] = ShortestPathTree(this, source, targets, dir);
            [distances, ix] = sort(distances, 'ascend');
            ordered_nodes = targets(ix);
        end
        
    end
    
    methods (Static)
        %% Verify Graph Connectivity
        % VERIFY_CONNECTIVITY verify if a graph is all connected. this method is
        % analogous to Width First Search. There, if we can construct a spanning tree
        % from any of the nodes in the graph, we think this graph is all connected.
        %
        %   [v, p_visit] = Graph.VerifyConnectivity(G, node_id)
        %
        % |G|: ajacent matrix of a graph.
        % |node_id|: a row vector, when this argument is specified, only check the
        % connectivity of this node set to other nodes.
        % |v|: true if the graph is all-connected.
        function [v, p_visit] = VerifyConnectivity(A, node_id)
            v = true;
            N = size(A,1);
            if N == 1
                error('error: Wrong Type of input argument G.');
            end
            if nargin == 1
                node_id = 1:N;
            end
            
            for n = node_id
                p_visit = zeros(1,N);   % nodes that have determined its neighbours if corresponding
                % element is 1, otherwise 0;
                t_visit = zeros(1,N);   % nodes that have been connected, but not fully determine its neighbours
                t_visit(n) = 1;
                i = find(t_visit,1);                % current visit node
                while isempty(i)==false
                    t_visit(i) = 0;
                    p_visit(i) = 1;
                    neighbors = find(A(i,:)~=0);
                    t = p_visit(neighbors)==0;      % new neighbours may be permanent nodes, which
                    % sould be eliminated from the find result
                    t_visit(neighbors(t)) = 1;
                    i = find(t_visit,1);            % current visit node
                end
                %% TODO
                % if only part of destination nodes should be ensured to connected to the srource
                % nodes, then p_visit only record these destination nodes.
                if ~isempty(find(p_visit==0,1))
                    v = false;
                    break;
                end
            end
        end
        %% random graph
        % the graph should be all connected. we can use width/depth first search to build
        % this random graph.
        function A = RandomAdjacent(N)
            A = zeros(N);
            degree_max = min(4, N-1);
            degree_min = 1;
            degree_ave = (degree_max+degree_min)/2;     % average degree is less than 3;
            neighbor_count = zeros(1,N);
            
            p_visit = zeros(1,N);       % nodes that have determined its neighbours if
            % corresponding element is 1, otherwise 0;
            t_visit = zeros(1,N);       % nodes that have been connected, but not fully determine
            t_visit(1) = 1;             % its neighbours
            
            while ~isempty(find(p_visit==0,1))
                i = find(t_visit,1);    % current visit node
                while isempty(i)==false
                    t_visit(i) = 0;
                    p_visit(i) = 1;
                    neighbor_num = randi([degree_min,degree_max]);
                    neighbors = find(A(i,:)~=0);
                    new_neighbors_num = neighbor_num - length(neighbors);
                    while new_neighbors_num > 0
                        potential_neighbor = ...
                            find((A(i,:)==0) & (neighbor_count < degree_ave) & (p_visit==0));
                        % find neighbours that have not been i's neighbour and degree less than
                        % average degree and have not been fully determined its neighbours;
                        if length(potential_neighbor) < new_neighbors_num
                            potential_neighbor = find((A(i,:)==0) & (neighbor_count < degree_max));
                            potential_neighbor(potential_neighbor==i)=[];
                        end
                        new_index = unique(randi(length(potential_neighbor),2*new_neighbors_num,1));
                        new_index = new_index(1:min(new_neighbors_num,length(new_index)));
                        new_neighbors = potential_neighbor(new_index);
                        A(i,new_neighbors) = 1;
                        neighbor_count(new_neighbors) = neighbor_count(new_neighbors)+1;
                        neighbor_count(i) = neighbor_count(i)+length(new_neighbors);
                        new_neighbors_num = new_neighbors_num - length(new_neighbors);
                        t = p_visit(new_neighbors)==0; % new neighbours may be permanent nodes
                        t_visit(new_neighbors(t)) = 1;
                    end
                    A(:,i) = A(i,:)';
                    i = find(t_visit,1);    % current visit node
                end
                i = find(p_visit==0,1);     % current visit node
                % link the node to the permanent node set
                if isempty(i)==false
                    t_visit(i) = 1;
                    potential_neighbor = find(p_visit);		% node i is not linking to any other nodes;
                    new_index = randi(length(potential_neighbor));
                    new_neighbors = potential_neighbor(new_index);
                    A(i,new_neighbors) = 1;
                    neighbor_count(new_neighbors) = neighbor_count(new_neighbors)+1;
                    neighbor_count(i) = 1;
                end
            end
            
            if ~isempty(find(A~=A',1))       % verification of symmetry
                error('error: the generated graph is not bi-directional symmetric.');
            end
            % Assertion: the graph should be all connected.
            if ~Graph.VerifyConnectivity(A)
                error('error: the generated graph is not all connected.');
            end
            fprintf('graph information:\n');
            fprintf('\t average degree: %f\n',mean(sum(A)));
            fprintf('\t min degree: %d\n', min(sum(A)));
        end
        
        function [sub_G, node_index] = SubGraph(A, src)
            % find the G's subgraph that contains node src.
            [v, p_visit] = Graph.VerifyConnectivity(A, src);
            if v
                sub_G = A;
            else
                node_index = find(p_visit==1);
                sub_G = A(node_index, node_index);
            end
        end
    end
    
end

%% See Also
%   GeoGraph < Graph
%   properties: node_names, location