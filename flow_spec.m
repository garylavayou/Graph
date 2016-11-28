function FlowTable = flow_spec(model, name, seed, graph, load_factor, flow_density)
%% Flow specifications
% function FlowTable = flow_spec(model='sample', name)
% function FlowTable = flow_spec(model='random', name, seed, graph)
%   name = 'flat', 'converge', 'fattree'
% function FlowTable = flow_spec(model='random', name, seed, graph, load_factor, flow_density)
%   name = '2-layer', '2-layer-cluster-hotspot', '2-layer-gateway-hotspot' 
if strcmp(model, 'sample')==true
    %% sample flow
    % 1 and 2 must be used corresponding to the graph mode (sample,1/2). the graph mode
    % (sample,1/2) can also use other traffic model(eg. random).
    if strcmp(name,'1')==true
        FlowTable = [
            1	4   2	1;
            1   5   1.5	2;
            2	6   1	3;
            3	4   0.5	4;
            3	5   1.2	5;
            4	6   2	6;
            4	1   3	7;
            5   6   1	8];
    elseif strcmp(name,'2')==true        % sample flow for G(15)
        FlowTable = [
            1  4   0.2	1;
            1  5   0.24	2;
            1  11  0.3	3;
            1  10  0.26	4;
            1  13  0.33	5;
            2  6   0.25	6;
            2  1   0.3	7;
            2  5   0.25	8;
            2  9   0.36	9;
            2  13  0.26	10;
            3  4   0.25	11;
            3  5   0.3	12;
            3  10  0.2	13;
            3  13  0.24	14;
            4  6   0.32	15;
            4  1   0.23	16;
            4  8   0.33	17;
            4  9   0.2	18;
            4  12  0.1	19;
            5  6   0.25	20;
            5  3   0.28	21;
            5  8   0.15	22;
            5  10  0.25	23;
            7  2   0.32	24;
            7  5   0.35	25;
            7  9   0.2	26;
            7  12  0.3	27;
            8  1   0.23	28;
            8  5   0.2	29;
            8  12  0.1	30;
            8  15  0.3	31;
            8  2   0.2	32;
            9  6   0.25	33;
            9  10  0.2	34;
            9  14  0.31	35;
            9  1   0.22	36;
            10 1   0.2	37;
            10 4   0.2	38;
            10 8   0.3	39;
            10 15  0.35	40;
            11 2   0.3	41;
            11 7   0.2	42;
            11 12  0.2	43;
            11 14  0.3	44;
            12 6   0.25	45;
            12 10  0.2	46;
            12 3   0.25	47;
            12 7   0.14	48;
            13 3   0.35	49;
            13 9   0.2	50;
            13 4   0.22	51;
            13 12  0.24	52;
            14 10  0.2	53;
            14 7   0.3	54;
            14 6   0.24	55;
            14 1   0.17	56;
            14 12  0.35	57;
            15 3   0.2	58;
            15 5   0.2	59;
            15 9   0.2	60;
            15 10  0.25	61];
    else
        error('Error: Invalid value of parameter name.');
    end
elseif strcmp(model,'random')==true
    demand_max = 10;            % tunable parameter
    demand_min = 1;
    rng(seed);
    
    if strcmp(name,'flat')==true
        warning('function removed.');
    elseif strcmp(name,'converge')==true
		warning('function removed.');
    elseif strcmpi(name, 'fattree')==true
		warning('function removed.');
    elseif strncmp(name,'2-layer',7)==true
        % Type-1 flow: from serving nodes to base stations.
        % Type-2 flow: between serving nodes.
        if nargin < 5 
            load_factor = 1;
        end
        
        num_core_nodes = graph.NumCoreNodes;
        if nargin < 6
            flow_density = num_core_nodes;
        end
        K = graph.Size*num_core_nodes;
        m = 6; v = 12;
        mu = log(m^2/sqrt(v+m^2));
        sigma = sqrt(log(v/m^2+1));
        
        FlowTable = zeros(K,4);
        FlowTable(:,4) = (1:K)';
        K = 0;
        % generate core network traffic (Type-2 Flow);
        if load_factor ~= 0
            for dest = 1:num_core_nodes
                src = unique(randi(num_core_nodes,min(flow_density,num_core_nodes),1));
                src(src==dest)=[];  % exclude the same src and dest.
                ls =length(src);
                FlowTable(K+1:K+ls,1)=src;
                FlowTable(K+1:K+ls,2)=dest;
                FlowTable(K+1:K+ls,3)=...
                    round(unifrnd(demand_min,demand_max,ls,1))*load_factor;
                K = K+ls;
            end
        end
        % generate cross-domain traffic (Type-1 Flow)
        if strcmp('2-layer-hotspot-gateway',name)==true
            for core_id = 1:num_core_nodes
                % generate access network traffic (Type-1 Flow);
                % find all node bind to the node(core_id)
                bind_acc_node = find(graph.AccessCoreMap==core_id)+num_core_nodes;
                bind_acc_node = reshape(bind_acc_node, 1, length(bind_acc_node));
                hotspot_ratio = lognrnd(mu,sigma);
                if hotspot_ratio > 18
                    hotspot_ratio = 18;
                elseif hotspot_ratio < 1
                    hotspot_ratio = 1;
                end
                for dest = bind_acc_node
                    % source of service only locate in core network.
                    src = unique(randi(num_core_nodes,min(flow_density,num_core_nodes),1));
                    ls =length(src);
                    FlowTable(K+1:K+ls,1)=src;
                    FlowTable(K+1:K+ls,2)=dest;
                    FlowTable(K+1:K+ls,3)=...
                        round(unifrnd(demand_min,demand_max,ls,1)*hotspot_ratio);
                    K = K+ls;
                end
            end
        elseif strcmp('2-layer-hotspot-cluster',name)==true || strcmp('2-layer',name)==true
            for dest = (num_core_nodes+1):graph.Size
                src = unique(randi(num_core_nodes,min(flow_density,num_core_nodes),1));
                src(src==dest)=[];  % exclude the same src and dest.
                ls =length(src);
                FlowTable(K+1:K+ls,1)=src;
                FlowTable(K+1:K+ls,2)=dest;
                if strcmp('2-layer-hotspot-cluster',name)==true
                    hotspot_ratio = lognrnd(mu,sigma);
                    if hotspot_ratio > 18
                        hotspot_ratio = 18;
                    elseif hotspot_ratio < 1
                        hotspot_ratio = 1;
                    end
                else
                    hotspot_ratio = 1;
                end
                FlowTable(K+1:K+ls,3) = ...
                    round(unifrnd(demand_min,demand_max,ls,1)*hotspot_ratio);
                K = K+ls;
            end
        else
            error('error: invalid flow model.');
        end
        FlowTable = FlowTable(1:K,:);
    end
else
    error('Error: Invalid flow model.');
end
end