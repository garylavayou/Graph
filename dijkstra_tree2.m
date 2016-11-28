function path_list = dijkstra_tree2(W,source,dest_set)
%DIJKSTRA_TREE2 generate a shortest path tree by dijkstra algorithm. this
%function is the extended version of DIJKSTRA_TREE, since one shortest path is
%find between one source and a group of destination.
% W: Directed connected graph;
% source: Source node;
% dest_set: the set of destination nodes organized as cell vector. Each cell
% element is a vector represents multiple destination.
% path_list: The tree is represented by a series of paths organised as a cell
% vector. Corresponding to the node in dest_set, each cell element is a path
% consists of a sequence of nodes from source to destination.
%NOTE: there are more efficient algorithms using fibonacci heap.
%% assertion
len_dest_set = length(dest_set);
N = size(W,1);		% the number of nodes in the graph
for i = 1:len_dest_set
    if sum(find(dest_set{i}==source))>0
        error('argument error: destination must be different from source.');
    end
    if sum(find(dest_set{i}>N))>0
        error('argument error: invalid destination node.');
    end
end
%% parameter declaration
dist = W(source,:);     	% distance to source node
b_permanent = zeros(1,N);	% permanent label for nodes
b_permanent(source) = 1;
prev = ones(1,N)*source;    % the recode of previous node on the shortest path

%% dijkstra algorithm
for i = 1:N-1
    temp_dist = dist;
    temp_dist(b_permanent==1) = inf;  
%     for j = 1:N
%         if b_permanent(j)
%             temp_dist(j) = inf;		
%         % else temp_dist(j) = dist(j);
%         end
%     end
    % find a nearest node to source
    % if a node is permanently labelled, do not consider it.
    % here, a connected graph is assumed, otherwise an assertion is needed to break the loop
    [~,index] = min(temp_dist);
    b_permanent(index) = 1;
    
    % update the distance function
    for k=1:N
        if dist(k) > dist(index)+W(index,k)
            dist(k) = dist(index)+W(index,k);
            prev(k) = index;
        end
    end 
end

%% back trace the shortest path
path_list = cell(1,length(dest_set));
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
    path_list(k) = {path(j:-1:1)};
    k = k + 1;
end
end