function capacity_list = max_capacity(W,source,dest_set)
%MAX_CAPACITY find the maximum capacity between source and each destination by dijkstra
%algorithm. One maximum capacity path is find from "source" to "destination" by
%DIJKSTRA_CAPACITY. 
% W: Directed connected graph, W(i,j) denote the edge capacity from node i to node j;
% source: Source node;
% dest_set: A vector represent the set of destination nodes. Each time, the
% maximum capacity path is retrived between source and one node from dest_set.
% capacity_list: Corresponding to the node in dest_set, each array element is the maximum
% capacity from source to destination.  
%% assertion
if sum(find(dest_set==source))>0
    error('argument error: destination must be different from source.');
end
N = size(W,1);		% the number of nodes in the graph
if max(dest_set)>N || source > N
    error('argument error: invalid node index.');
end  

%% Dijkstra Algorithm
[~, capacity] = dijkstra_capacity(W,source);
capacity_list = capacity(dest_set);
end