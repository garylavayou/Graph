function [path_set, new_path_index] = path_set_merge(path_set, path_list)
%path_set_merge merge PATH_LIST into PATH_SET.
% path_set: (N*N) cell array, in which each element {i,j} is also a cell
%                   array, which reprensent the list of path from i to j;
% path_list: the list of path to be merged.
new_path_index = zeros(length(path_list),3);
for i = 1:length(path_list)
    path = path_list{i};
    src = path(1);
    dest = path(length(path));
    index = path_index(path_set,path);
    if isempty(index)   % no path exists
        path_set{src,dest} = struct('path',path,'bandwidth',0);
        new_path_index(i,:)=[src dest 1];
    elseif sum(index)==0  % this path is different from existing path
        path_num = length(path_set{src,dest});
        path_set{src,dest}(path_num+1)=struct('path',path,'bandwidth',0);
        new_path_index(i,:)=[src dest path_num+1];
    else % this path already exist   
        new_path_index(i,:)=index;
    end
end
end