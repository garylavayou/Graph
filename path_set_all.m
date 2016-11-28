%% path_set_all: 
% *IT IS IMPOSSIBLE TO ATTAIN ALL PATHS EVEN IN A MEDIUM SIZE GRAPH* 
%% Function declaration
% function [path_set,path_num] = path_set_all(W)
function [path_set,path_num] = path_set_all(W)  %, st_pair
N = size(W,1);
adj_list = make_adj_list(W);
% st_matrix = zeros(N);
% for i=1:length(st_pair)
%     st_matrix(st_pair(i,1),st_pair(i,2))=1;
% end
%% start visit
hwait=waitbar(0,'Please wait >>>>>>>>');
path_set = cell(N);
progress_all = N^2;
progress = 0;
for source = 1:N   % visit all path from source node i(i=1,...,N)
    b_visit = zeros(1,N);
    top = 0;        % counter for stack top
    stack = zeros(1,N);
    v = source; b_visit(source) = 1;
    top = top + 1; stack(top) = v;
    head = 1;    % flag to visit adj_list
    p = 0; neighbors = 0;
    while top>0
        if head == 1
            neighbors = adj_list{v};
            p = 1;
            head = 0;
        else
            p = p + 1;
        end
        if p<=length(neighbors)
            if b_visit(neighbors(p))==0
                b_visit(neighbors(p)) = 1;
                top = top + 1; stack(top) = neighbors(p);
                v = stack(top);
                head = 1;
            end
        else
            % output a path
            if isempty(path_set{stack(1),stack(top)})
                path_set{stack(1),stack(top)} = {stack(1:top)};   % path_set(i,i) is invlalid
                progress = progress + 1;
                waitbar(progress/progress_all,hwait,'Almost done');
            else
                len = length(path_set{stack(1),stack(top)});
                path_set{stack(1),stack(top)} = {path_set{stack(1),stack(top)}{1:len} stack(1:top)};
            end
            b_visit(stack(top)) = 0; top = top - 1;
            if top > 0  % 定位下一个未遍历的节点
                neighbors = adj_list{stack(top)};
                p = 1;
                while neighbors(p) ~= v
                    p = p + 1;
                end
                v = stack(top);
                head = 0;
            end
        end
    end
end
close(hwait);
path_num = 0;
for i = 1:N
    for j = 1:N
       path_num = path_num + length(path_set{i,j}); 
    end
end
end

function adj_list = make_adj_list(W)
N = size(W,1);
adj_list = cell(1,N);
for i =1:N;
    node_list = zeros(1,N);
    k = 0;
    for j =1:N
        if W(i,j)~=0
            k = k+1;
            node_list(k) = j;
        end
    end
    adj_list(i)={node_list(1:k)};
end
end