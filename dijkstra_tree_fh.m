function path_list = dijkstra_tree_fh(fh,Gl,source,dest_set, dir)
%DIJKSTRA_TREE generate a shortest path tree, which consists of shortest paths from
%"source" to each "destination", by dijkstra algorithm. It also can calculate the reverse
%shortest path tree, whose paths are from destinations to source.
% W: Directed all-connected graph; W(i,j) denote the distance from node i to node j;
% source: Source node;
% dest_set: A vector represent the set of destination nodes. the shortest path
% is retrieved between source and one node from dest_set. 
% dir: 0 - shortest paths from source to all other nodes(default);
%      1 - shortest paths from all other nodes to source;
% path_list: The tree is represented by a series of paths organised as a cell
% vector. Corresponding to the node in dest_set, each cell element is a path
% consists of a sequence of nodes from source to destination. 
% NOTE: this implementation use dijkstra algorithm based on Fibonacci heap.
%% assertion
N = size(Gl,1);		% the number of nodes in the graph
switch nargin
    case 3
        dest_set = (1:N)';
        dest_set(source)=[];
        dir = 0;
    case 4
        dir = 0;
end
if sum(find(dest_set==source))>0
    error('argument error: destination must be different from source.');
end
if max(dest_set)>N || source > N
    error('argument error: invalid node index.');
end  
if dir
    Gl = Gl';                 % reverse the direction of edges in the graph
end

%% Dijkstra Algorithm
dist = Gl(single_src,:);     	% distance from source node to other joint nodes
prev = ones(N,1)*single_src;    % the recode of previous node on the shortest path

for v = 1:N
    fh.insert(v, dist(v));      % dist(src)=0, dist(src,neighbor)=l(e)
                                % dist(src,non-neighbor)=inf
end

while ~isempty(fh.min_node)
    closest_node = fh.delete_min();
    u = closest_node.item;
    for v = find(Gl(u,:)<Inf)
        tent_dist = dist(u) + Gl(u,v);
        if tent_dist < dist(v)
    		fh.decrease_key(v, tent_dist);
			dist(v) = tent_dist;
			prev(v) = u;
        end
    end
end


%% (back) trace the shortest path
path_list = cell(length(dest_set),1);
k = 1;
for i = dest_set'   % i~=source
    path = zeros(1,N);
    path(1) = i;
    j = 1;
    while path(j) ~= source
        t = prev(path(j));
        j = j+1;
        path(j) = t;
    end
    if dir
        %  Since the direction is dest to src, backtrace is not required.
        path_list(k) = {path(1:j)};
    else
        path_list(k) = {path(j:-1:1)};
    end
    k = k + 1;
end
end