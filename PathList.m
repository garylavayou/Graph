%% Path List
% Candidate paths of a flow between two nodes.
classdef PathList < matlab.mixin.Copyable
    
    properties
        Source;
        Destination;
        Width;
        Length;
        Latency;
		end
		
		properties (Access = protected)				
				path_list;
		end
		
		properties (Constant, Access = private)
			storage_class = ?Path;
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
					this.path_list = ListArray('Path');
					if nargin == 0
						return;
					elseif isa(path_list, 'PathList')
						this = path_list.copy;
					else
						this.Add(path_list);
					end
        end
        
        function delete(this)
            delete(this.path_list);
        end
		end
		    
    methods
        function n = get.Width(this)
            n = this.path_list.Length;
        end
        
        function l = get.Length(this)
					l = max(this.path_list{'Length'});
				end
        
        function l = get.Latency(this)
					l = max(this.path_list{'latency'});
        end
        
        function s = get.Source(this)
					if this.Width == 0
						s = [];
					else
            s = this.path_list{1}.Source;
					end
        end
        
        function t = get.Destination(this)
					if this.Width == 0
						t = [];
					else
						t = this.path_list{1}.Destination;
					end
				end
				
				function n = numArgumentsFromSubscript(this,s,indexingContext) %#ok<INUSL>
					%%%
					% see also <ListArray>.
					switch indexingContext
						case matlab.mixin.util.IndexingContext.Statement
							n = 0;
						case matlab.mixin.util.IndexingContext.Expression
							n = 1; % nargout for indexed reference used as function argument
						case matlab.mixin.util.IndexingContext.Assignment
							n = 1; % nargin for indexed assignment
					end
				end
				
				%% Deprecated
				% cannot distinguish the {} and () operator.
				%
				% 				function ind = end(this,~,~)
				%             ind = this.path_list.Length;
				%         end
				
				%%
				% see also <ListArray>.
				function varargout = subsref(list, s)
					switch s(1).type
						case '.'
							assertpermission('PathList', s(1).subs, 'get');	% check the access permission
							dims = size(list);
							elements = cell(dims);
							for i = 1:numel(list)
								elements(i) = {builtin('subsref',list(i),s)};
							end
							b_concat = assertcat(elements, true);
							varargout{1} = tryconcat(elements, dims, b_concat);
						case '()'
							sub_list = list(s(1).subs{:});
							if length(s) == 1
								varargout{1} = sub_list;
								return;
							elseif ~isequal(s(2).type, '.')
								error('error: only support obj(subs).name operation.');
							end
							assertpermission('PathList', s(2).subs, 'get');	% check the access permission
							dims = size(sub_list);
							elements = cell(dims);
							for i = 1:numel(sub_list)
								elements(i) = {builtin('subsref', sub_list(i), s(2:end))};
							end
							b_concat = assertcat(elements, true);
							varargout{1} = tryconcat(elements, dims, b_concat);
						case '{}'
							if ~isscalar(list)
								error('error: ''{}'' operation only supported for scalar <PriorityQueue>.');
							end
							if length(s) == 1
								if length(s(1).subs) == 1 && ischar(s(1).subs)
									assertpermission('Path', s(1).subs, 'get');
								elseif length(s(1).subs) == 2 && ischar(s(1).subs{2})
									assertpermission('Path', s(1).subs{2}, 'get');
								end
							else
								assertpermission('Path', s(2).subs, 'get');
							end
							varargout = {builtin('subsref', list.path_list, s)};
						otherwise
							error('error: operation %s is not supported.', s(1).type);
					end
				end
				
				function list = subsasgn(list, s, v)
					switch s(1).type
						case '.'
							sub_list = list;
						case '()'  % v is cell array
							if isempty(list) || length(s) == 1
								list = builtin('subsasgn', list, s, v);
								return;
							end
							sub_list = list(s(1).subs);
							s = s(2:end);
						case '{}'
							if ~isscalar(list)
								error('error: ''{}'' operation only supported for scalar <PriorityQueue>.');
							end
							if length(s) == 1
								if length(s(1).subs) == 1 && ischar(s(1).subs)
									assertpermission('Path', s(1).subs, 'set');
								elseif length(s(1).subs) == 2 && ischar(s(1).subs{2})
									assertpermission('Path', s(1).subs{2}, 'set');
								end
							else
								assertpermission('Path', s(2).subs, 'set');
							end
							list.path_list = builtin('subsasgn', list.path_list, s, v);		
							return;
						otherwise
							error('error: operation %s is not supported.', s(1).type);
					end
					
					assertpermission('PathList', s(1).subs, 'set');
					idx = 1:numel(sub_list);
					if ischar(v)
						v = {v};
					end
					lenv = numel(v);
					if lenv~=1 && lenv~=length(idx)
						error('error: the number of required value is inconsistent with the supply.');
					end
					val = cell(size(idx));
					for i = length(idx)
						if isempty(v)
							val{i} = [];
						else
							if iscell(v) || isstring(v)
								val{i} = v{min(lenv, i)};
							else
								val{i} = v(min(lenv, i));
							end
						end
					end
					for i = 1:numel(idx)
						sub_list(i) = builtin('subsasgn', sub_list(i), s, val{i});  % list is handle, the value change takes effect
					end
				end
		end
    
		methods
			function p = Add(this, path_list)
				p = Path.empty();
				if isnumeric(path_list)
					p = this.path_list.Add(Path(path_list));
				elseif isa(path_list, 'Path')
					p = this.path_list.Add(path_list);
				elseif iscell(path_list)
					for i = 1:length(path_list)
						if i> 1 && ~PathList.assert_common_path(path_list{i}, path_list{i-1})
							error('[%s] error: Paths have different source/destination.', calledby(0));
						end
						if isnumeric(path_list{i})
							p = this.path_list.Add(Path(path_list{i}));
						else
							p = this.path_list.Add(path_list{i});
						end
					end
				else
					error('[%s] error: invalid arguments.', calledby(0));
				end
			end
		end
		
		methods (Access = protected)
			function obj = copyElement(this)
				obj = copyElement@matlab.mixin.Copyable(this);
				obj.path_list = this.path_list.copy();
			end
			
			function idx = assertindex(this, indices)
				% See also <ListArray>
				if isempty(indices)
					idx = 0;
				elseif ischar(indices)
					if isequal(indices, ':')
						idx = 1:this.path_list.Length;
					else
						error('error: invalid index')
					end
				else
					if islogical(indices)
						indices = find(indices);
					end
					if max(indices) > this.path_list.Length
						error('error: Index out of bound.');
					end
					if min(indices) <= 0
						error('error: negative index.');
					end
					idx = indices;
				end
			end
		end
		
		methods (Access = protected, Static)
			function tf = assert_common_paths(p1, p2)
				if isempty(p1) || isempty(p2)
					tf = true;
					return;
				end
				if isnumeric(p1)
					p1 = Path(p1);
				end
				if isnumeric(p2)
					p2 = Path(p2);
				end
				
				tf = p1.Source == p2.Source && p1.Destination == p2.Destination;
			end
		end
end

