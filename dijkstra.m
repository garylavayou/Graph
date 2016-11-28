function [prev, dist] = dijkstra(W,source,dest_set, dir)
%DIJKSTRA generate a shortest path tree, which consists of shortest paths from
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
% NOTE: there are more efficient algorithms using Fibonacci heap.
%% assertion
N = size(W,1);		% the number of nodes in the graph
switch nargin
    case 2
        dest_set = (1:N)';
        dest_set(source)=[];
        dir = 0;
    case 3
        dir = 0;
end
if sum(find(dest_set==source))>0
    error('argument error: destination must be different from source.');
end
if max(dest_set)>N || source > N
    error('argument error: invalid node index.');
end  
if dir
    W = W';                 % reverse the direction of edges in the graph
end

%% Dijkstra Algorithm
dist = W(source,:)';     	% distance from source node to other joint nodes
                            % equal to W(:,source) if dir=1
b_permanent = zeros(1,N);	% permanent label for nodes
b_permanent(source) = 1;
prev = ones(N,1)*source;    % the recode of previous node on the shortest path
% count = length(dest_set);

% O(N) = 5N^2 +14N
for i = 1:N-1
    temp_dist = dist;
    temp_dist(b_permanent==1) = inf;  
    % find a nearest node to source
    % if a node is permanently labelled, do not consider it.
    % here, a connected graph is assumed, otherwise an assertion is needed to break the loop
    [min_dist,index] = min(temp_dist);
    if min_dist == inf
        error('Error: the graph is not all connected.');
    end
    b_permanent(index) = 1;
    % update the distance function
    for k = find(W(index,:)<Inf)
        % Permanent nodes' distance and previous node will not change any more,
        % so they are excluded from the distance update procedure.
        t = dist(index)+W(index,k);         % W(u,u) = 0;
        if dist(k) > t
            dist(k) = t;
            prev(k) = index;
            % For the node added to permanent set just above, its previous node
            % has been determined in the previous iteration or in the
            % initialization phase(source node is its previous node).
        end
    end 
    % If the dest_set is always far less than the node set, this assertion is necessary.
%     if ~isempty(find(dest_set==index,1))
%         count = count - 1;
%         if count == 0
%             break;
%         end
%     end
end
end