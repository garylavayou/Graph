classdef TwoLayerGraph < Graph
    %TWOLAYERGRAPH the graph comprises two layers of random graph.
    %   The connection between the two layer is randomly decided.
    
    properties
        upper_size;         % the size of the upper layer graph
        bottom_size;        % the size of the bottom latyer graph
        
    end
    
    methods
        %% Constructor
        %   Build(network)
        %   Build(layer_size, seed, link_capacity)
        %
        % |layer_size|: the graph size of each layer. If layer_size is scalar, the two
        %     layer have the same size.
        %
        % |seed|: seed of random number generator for creating graph.
        %
        % |link_capacity|: the link capacity of each layer. If this argument is not
        %     specified, the link capcity is set to 1. If |link_capcity| is a scalar,
        %     the link capacity of the two layer is the same.
        function this = TwoLayerGraph(arg1, seed, link_capacity)
            if nargin == 1
                network = arg1;
                this.bottom_size = network.BottomSize;
                this.upper_size = network.UpperSize;
                graph_size = this.upper_size + this.bottom_size;
                this.Adjacent = zeros(graph_size, graph_size);
                for i = 1:network.IndirectedLinkNumber
                    this.Adjacent(network.link(i,1), network.link(i,2)) = ...
                        network.link(i,4);
                    this.Adjacent(network.link(i,2), network.link(i,1)) = ...
                        network.link(i,4);
                end
            else
                layer_size = arg1;
                if nargin >= 2
                    rng(seed);
                else
                    rng(20160519);
                end
                if nargin >= 3
                    if numel(link_capacity) == 1
                        link_capacity =  ones(2,1) * link_capacity;
                    end
                else
                    link_capacity = [1;1];
                end
                
                if numel(layer_size) == 1
                    layer_size = ones(2,1)*layer_size;
                end
                this.upper_size = layer_size(1);
                this.bottom_size = layer_size(2);
                graph_size = sum(layer_size);
                this.Adjacent = zeros(graph_size);
                
                % generate upper layer
                node_index = 1:this.upper_size;
                this.Adjacent(node_index, node_index) = ...
                    Graph.RandomAdjacent(this.upper_size) * link_capacity(1);
                % generate bottom layer
                node_index = this.upper_size + (1:this.bottom_size);
                this.Adjacent(node_index, node_index) = ...
                    Graph.RandomAdjacent(this.bottom_size) * link_capacity(2);
                % TODO: determine the connection between two layer
            end
            
            this.Renew;      % Update the properties in the object.
        end
    end  
end