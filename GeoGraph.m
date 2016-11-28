classdef GeoGraph < Graph
    
    properties
        node_names;
        locations;
    end
    
    methods        
        function Copy(this, graph)
            Copy@Graph(graph);
            this.node_names = graph.node_names;
            this.locations = graph.locations;
            this.distance = graph.distance;
        end
    end
    
    methods(Access=protected)
        function Renew(this)
            Renew@Graph;
            this.distance = 1./this.Adjacent;
            for i = 1:size(this.Adjacent,1)
                this.Distance(i,i) = 0;
            end
        end
    end
end

