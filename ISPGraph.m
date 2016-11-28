classdef ISPGraph < Graph
    %ISPGRAPH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        node_id;
        node_name;
        node_pos;
        num_origin_edge;
		
		NodeWeight;
    end
    
    methods
        function this = ISPGraph(filename)
			if nargin == 0
				return;
			end
		
            GraphDoc = xmlread(filename);
            node_list = GraphDoc.getElementsByTagName('node');
            n = node_list.getLength;
            this.Size = 0;
            this.node_id = zeros(n, 1);
            this.node_name = cell(n, 1);
            this.node_pos = zeros(n, 2);
            for k=0:n-1
                node = node_list.item(k);
                node_data_list = node.getElementsByTagName('data');
                internal = 1;
                for j = 0:node_data_list.getLength-1
                    data = node_data_list.item(j);
                    key = data.getAttribute('key');
                    if strcmp('Internal',key)==true
                        internal = str2double(data.getTextContent);
                        break;
                    end
                end
                if internal == 1
                    this.Size = this.Size + 1;
                    i = this.Size;
                    this.node_id(i) = str2double(node.getAttribute('id'));
                    for j = 0:node_data_list.getLength-1
                        data = node_data_list.item(j);
                        key = data.getAttribute('key');
                        if strcmp('label',key)==true
                            this.node_name{i} = char(data.getTextContent);
                            continue;
                        end
                        if strcmp('Longitude',key)==true
                            this.node_pos(i,1) = str2double(data.getTextContent);
                            continue;
                        end
                        if strcmp('Latitude',key)==true
                            this.node_pos(i,2) = str2double(data.getTextContent);
                            continue;
                        end
                    end
                end
            end
            this.node_id(this.Size+1:n) = [];
            this.node_name(this.Size+1:n) = [];
            this.node_pos(this.Size+1:n,:) = [];
            
            edge_list = GraphDoc.getElementsByTagName('edge');
            n = edge_list.getLength;
            this.num_origin_edge = 0;
            this.Adjacent = zeros(this.Size);
            for k = 0:n-1
                edge = edge_list.item(k);
                source = int32(str2double(edge.getAttribute('source')));
                dest = int32(str2double(edge.getAttribute('target')));
                s_index = find(source==this.node_id,1);
                d_index = find(dest==this.node_id,1);
                if ~isempty(s_index)&& ~isempty(d_index)
                    edge_data_list = edge.getElementsByTagName('data');
                    for j = 0:edge_data_list.getLength-1
                        data = edge_data_list.item(j);
                        key = data.getAttribute('key');
                        if strcmp('LinkSpeedRaw',key)==true
                            % there may exist multiple links between two nodes, so aggregating it as
                            % one link.
                            this.Adjacent(s_index,d_index) = ...
                                this.Adjacent(s_index,d_index) + str2double(data.getTextContent);
                            continue;
                        end
                    end
                    this.num_origin_edge = this.num_origin_edge + 1;
                end
            end
            gh = GraphDoc.getElementsByTagName('graph');
            direct = gh.item(0).getAttribute('edgedefault');
            if strcmp(direct, 'undirected')==true
                this.Adjacent = this.Adjacent+this.Adjacent';
            end
            
			this.Renew;      % Update the properties in the object.
        end 
    end
    methods(Access=protected)
        function Renew(this)
            Renew@Graph;
            this.NodeWeight = sum(this.Adjacent);
			this.NodeWeight = this.NodeWeight'/min(this.NodeWeight);
        end
    end    
end

