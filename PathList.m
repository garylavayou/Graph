%% Path List
% Candidate paths of a flow between two nodes.

classdef PathList < handle
    
    properties
        Source;
        Destination;
        Width;
        Length;
        Latency;
        
        paths;
    end
    
    methods
        %% Constructor
        % Input Aregument
        %
        % |path_list|: a list of path instance of Path class.
        %     Note: Path is a handle class, so the constructor must recreate a copy of the
        %     input argument.
        %
        function this = PathList(path_list)
            if isa(path_list, 'PathList')
                this.copy(path_list);
            else
                lp = length(path_list);
                this.paths = cell(lp,1);
                for i= 1:lp
                    this.paths{i} = Path(path_list{i});
                end
            end
        end
        
        
    end
    
    methods
        function n = get.Width(this)
            n = length(this.paths);
        end
        
        function l = get.Length(this)
            l = 0;
            for i = 1:this.Width
                l = max(l,this.paths{i}.Length);
            end
        end
        
        function l = get.Latency(this)
            l = 0;
            for i = 1:this.Width
                l = max(l,this.paths{i}.latency);
            end
        end
        
        function s = get.Source(this)
            s = this.paths{1}.Destination;
        end
        
        function t = get.Destination(this)
            t = this.paths{1}.Source;
        end
    end
    
    methods (Access = private)
        function copy(this, path_list)
            lp = path_list.Width;
            this.paths = cell(lp,1);
            for i= 1:lp
                this.paths{i} = Path(path_list.paths{i});
            end
        end
    end
end

