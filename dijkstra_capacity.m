function [prev, capacity] = dijkstra_capacity( W,source, dir )
%DIJKSTRA_CAPACITY Calculate the maximum capacity paths from source to all other nodes or
%the reverse direction path. This function is distinguished from DIJKSTRA by the adjoint
%matrix whose elements stand for capacity of edge. 
% W: Directed connected graph, W(i,j) denote the edge capacity from node i to node j;
% source: Source node;
% dir: 0 - objective paths from source to all other nodes(default);
%      1 - objective paths from all other nodes to source;
% prev: the output is a vector which records a tree structure, prev(i) represent node i's
% parent node in the tree. Using this structure, the path from or to source can be
% traced(dir=1) or backtraced(dir=0); 
% NOTE: there are more efficient algorithms using fibonacci heap.

%% parameter declaration
if nargin == 2
    dir = 0;
end
if dir
    W = W';                 % reverse the direction of edges in the graph
end
N = size(W,1);              % the number of nodes in the graph
capacity = W(source,:)';    % capacity from source node to other adjoint nodes
                            % equal to W(:,source) if dir=1
b_permanent = zeros(N,1);	% permanent label for nodes
b_permanent(source) = 1;
prev = ones(N,1)*source;    % the recode of previous node on the shortest path

%% dijkstra algorithm
for i = 1:N-1
    temp_cap = capacity;
    temp_cap(b_permanent==1) = 0;
    % Find a temporary node which has the maximum capacity from/to the source
    % If a node is permanently labelled, do not consider it.
    % Here, a connected graph is assumed, otherwise an assertion is needed to break the
    % loop.
    [~,index] = max(temp_cap);
    b_permanent(index) = 1;
    
    % update the capacity function: c(p) = min(c(p1),c(p2)), p=p1+p2;
    for k=(find(b_permanent==0))'
        % Permanent nodes' capacity and previous node will not change any more,
        % so they are excluded from the capacity update procedure.
        cap_index_k = min(capacity(index),W(index,k));
        if capacity(k) < cap_index_k
            capacity(k) = cap_index_k;
            prev(k) = index;
            % For the node added to permenent set just above, its previous node
            % has been determinted in the previous iteration or in the
            % intilization phase(source node is its previous node).
        end
    end 
end

end

