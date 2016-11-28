%% Sample Graph
%   SampleGraph < Graph
%
% Smaple graphes for testing algorithm
%%
classdef SampleGraph < Graph    
    properties
    end
    methods
        function this = SampleGraph(name)
            if nargin == 0
                return;
            end
            % adjacent matrix with capacity
            if strcmpi(name,'sample1')==true
                this.Adjacent = [   % source:
                    0  7  9  0 0 14;
                    7  0 10 15 0  0;
                    9 10  0 11 0  2;
                    0 15 11  0 6  0;
                    0  0  0  6 0  9;
                    14 0  2  0 9  0];
            elseif strcmpi(name,'sample2')==true
                this.Adjacent = [   % source: traffic engineering in software defined network
                    0 1 1 1 0 0 0 0 0 0 0 0 0 0 0;    %1
                    1 0 1 0 1 0 0 0 0 0 1 0 0 0 0;    %2
                    1 1 0 1 0 1 1 0 0 0 0 0 0 0 0;    %3
                    1 0 1 0 0 0 0 1 1 0 0 0 0 0 0;    %4
                    0 1 0 0 0 0 0 0 0 0 0 1 0 0 0;    %5
                    0 0 1 0 0 0 1 0 0 1 1 0 0 0 0;    %6
                    0 0 1 0 0 1 0 0 1 1 0 0 0 0 0;    %7
                    0 0 0 1 0 0 0 0 1 0 0 0 0 0 0;    %8
                    0 0 0 1 0 0 1 1 0 1 0 0 0 0 1;    %9
                    0 0 0 0 0 1 1 0 1 0 1 0 1 1 0;    %10
                    0 1 0 0 0 1 0 0 0 1 0 1 1 0 0;    %11
                    0 0 0 0 1 0 0 0 0 0 1 0 1 0 0;    %12
                    0 0 0 0 0 0 0 0 0 1 1 1 0 1 0;    %13
                    0 0 0 0 0 0 0 0 0 1 0 0 1 0 1;    %14
                    0 0 0 0 0 0 0 0 1 0 0 0 0 1 0];   %15
            elseif strcmpi(name, 'sample3') == true
                this.Adjacent = [  % source: .\Joint TE\topology\random_network_10
                    0 1 0 0 1 0 1 1 0 0;    % n(1) = 2,5,7,8
                    1 0 0 1 0 0 1 0 0 1;    % n(2) = 1,4,7,10
                    0 0 0 1 0 1 0 0 0 0;    % n(3) = 4,6
                    0 1 1 0 1 0 0 0 0 0;    % n(4) = 2,3,5
                    1 0 0 1 0 0 1 1 0 0;    % n(5) = 1,4,7,8
                    0 0 1 0 0 0 0 0 0 1;    % n(6) = 3,10
                    1 1 0 0 1 0 0 0 0 0;    % n(7) = 1,2,5
                    1 0 0 0 1 0 0 0 1 1;    % n(8) = 1,5,9,10
                    0 0 0 0 0 0 0 1 0 1;    % n(9) = 8,10
                    0 1 0 0 0 1 0 1 1 0]*4;   % n(10)= 2,6,8,9
            else
                this.Adjacent = SampleGraph.LoadBackboneAdjacent(name);
            end
            
            this.Renew;     % Update the properties in the object.
        end
    end
    methods (Static)
        %% Method: Load backbone topology
        %   GA = LoadBackboneAdjacent(backbone_name)
        function GA = LoadBackboneAdjacent(backbone_name)
            % SmallNet
            if strcmpi('smallnet', backbone_name)==true
                net = [
                    0,1,0;
                    0,5,0;
                    0,6,0;
                    1,2,0;
                    1,6,0;
                    1,7,0;
                    2,3,0;
                    2,7,0;
                    2,8,0;
                    3,4,0;
                    3,8,0;
                    4,5,0;
                    4,8,0;
                    4,9,0;
                    5,6,0;
                    5,9,0;
                    6,7,0;
                    6,8,0;
                    6,9,0;
                    7,8,0;
                    8,9,0]+1;
            elseif strcmpi('ARPA2', backbone_name)==true
                net = [
                    0, 1,0;
                    0, 3,0;
                    0, 7,0;
                    1, 2,0;
                    2, 5,0;
                    3, 4,0;
                    4, 5,0;
                    5, 6,0;
                    5,14,0;
                    6, 7,0;
                    7, 8,0;
                    7,12,0;
                    8, 9,0;
                    9,10,0;
                    10,11,0;
                    10,16,0;
                    11,13,0;
                    12,13,0;
                    13,15,0;
                    14,15,0;
                    15,18,0;
                    16,17,0;
                    17,20,0;
                    18,19,0;
                    19,20,0]+1;
            elseif strcmpi('bellcore', backbone_name)==true
                net = [
                    0, 1,0;
                    0, 8,0;
                    0, 9,0;
                    1, 2,0;
                    1, 7,0;
                    1, 9,0;
                    1,10,0;
                    1,12,0;
                    2, 3,0;
                    2, 5,0;
                    2,12,0;
                    3, 4,0;
                    3,12,0;
                    4, 5,0;
                    4,14,0;
                    5, 6,0;
                    5,11,0;
                    5,13,0;
                    5,14,0;
                    6, 7,0;
                    6,11,0;
                    7, 8,0;
                    7,10,0;
                    7,11,0;
                    8, 9,0;
                    8,10,0;
                    11,12,0;
                    11,13,0]+1;
            elseif strcmpi('nsfnet', backbone_name)==true
                net = [
                    0, 1,0;
                    0, 2,0;
                    0, 3,0;
                    1, 2,0;
                    1, 7,0;
                    2, 5,0;
                    3, 4,0;
                    3, 9,0;
                    4, 5,0;
                    4, 6,0;
                    5,10,0;
                    5,11,0;
                    6, 7,0;
                    7, 8,0;
                    8,10,0;
                    8,12,0;
                    8,13,0;
                    9,12,0;
                    9,13,0;
                    11,12,0;
                    11,13,0;
                    ]+1;
            else
                error('Invalid backbone name.');
            end
            % backbone traffic is proportional to the destination nodes weights
            % since the important nodes will provide more service to other nodes.
            
            %% TODO data center topology

            %% TODO artificial core network topology
            
            num_node = max(max(net(:,1:2)));
            GA = zeros(num_node);
            for i = 1:size(net,1)
                GA(net(i,1),net(i,2)) = net(i,3);
            end
            GA = GA+GA';
        end
    end
    
end

