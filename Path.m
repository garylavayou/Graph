%% Path
% A path has a node list and associated with parameters, including, bandwidth, latency,
% identifier, etc.
%%
classdef Path < matlab.mixin.Copyable
    
    properties
        node_list;
        %% Bandwidth, Latency and Identifier
        % Bandwidth, latency and identifier are not the intrinsic properties of a path. So
        % these properties is assigned when the path is build from a graph and can be
        % altered thereafter.
        bandwidth;
        latency;
        id;
        
        Length;
        Source;
        Destination;
        TailLink;
        HeadLink;
    end
    
    methods

        function this = Path(node_list, bandwidth, latency)
            if nargin >= 1
                if isa(node_list, 'Path')   % node_list is an instance of Path class.
                    this = node_list.copy;
                elseif isnumeric(node_list)
                    this.node_list = node_list;
                else
                    error('error: invalid input argument (node_list).');
                end
            end
            if nargin >=2 && ~isempty(bandwidth)
                this.bandwidth = bandwidth;
            end
            if nargin >= 3 && ~isempty(latency)
                this.latency = latency;
            end
        end
        

        function l = get.Length(this)
            l = length(this.node_list);
        end
        function h = get.Source(this)
            h = this.node_list(1);
        end
        function t = get.Destination(this)
            t = this.node_list(end);
        end
        function e = get.TailLink(this)
            e = [this.node_list(end-1) this.node_list(end)];
        end
        function e = get.HeadLink(this)
            e = [this.node_list(1) this.node_list(2)];
        end
        
        function set.node_list(this, node_list)
            if isrow(node_list)
                this.node_list = node_list;
            else
                this.node_list = node_list';
            end
        end
        
        
        function n = Node(this, i)
            n = this.node_list(i);
        end
        function l = Link(this, i)
            l = [this.node_list(i) this.node_list(i+1)];
        end
        function b = EqualPath(this, path)
            if isequal(path.node_list, this.node_list)
                b = true;
            else
                b = false;
            end
        end
        
        function MultiplyBandwidth(this, factor)
            this.bandwidth = this.bandwidth*factor;
        end
        
        
    end
end

%% Constructor
% Build path with a list of node:
%
%      Path(node_list, bandwidth, latency)
%
% Build path with existing Path instance:
%
%      Path(path)
%

%% Property
% * *Length* |get|
% * *Source* |get|
% * *Destination* |get|
% * *TailLink* |get|
% * *HeadLink* |get|
% * *node_list* |set|

%% Methods
% * *Node*
%
%      n = Node(path, i)
%
% * *Link*
%
%      l = Link(path, i)
%
% * *EqualPath*
%
%      b = EqualPath(path1, path2)
