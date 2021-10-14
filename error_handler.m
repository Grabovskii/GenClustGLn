classdef error_handler
   
    properties
        n
        
        build_ghf
        build_c
        build_p
        bracket
        build_birat
    end
    
    methods
        function obj = error_handler(varargin)
            len = length(varargin);
            for i = 1:len
               if strcmpi(varargin{i},'n')
                   obj.n = varargin{i+1};
               elseif strcmpi(varargin{i},'build_ghf')
                   obj.build_ghf = varargin{i+1};
               elseif strcmpi(varargin{i},'build_c')
                   obj.build_c = varargin{i+1};
               elseif strcmpi(varargin{i},'build_p')
                   obj.build_p = varargin{i+1};
               elseif strcmpi(varargin{i},'bracket')
                   obj.bracket = varargin{i+1}; 
               elseif strcmpi(varargin{i},'build_birat')
                   obj.build_birat = varargin{i+1}; 
               end
            end
        end
        function output = bracket_check(lets,inds)
            %I assume that lets and inds are already converted
            %to avoid doing it twice
            output = cell(0); %if everything's okay, return empty output
            if length(lets) < 2 || length(inds) < 2
                output = {'The functions are not specified. Error.'};
                return;
            end
            
        end
        
        function b = check_existence(obj,l,i,varargin)
           %only checks whether a variable of type l has been constructed 
           %it won't check the indices, except for the lower row of phi
           if ~isempty(varargin)
               j = varargin{1};
           else
               j = cell(0);
           end

             if contains('ghfxy',l,'IgnoreCase',true) && obj.build_ghf
                 b = 1;
                 return;
             elseif strncmpi(l,'p',1) 
                 if ~(isempty(i)||isempty(j))
                    if ((i + j) == obj.n && obj.build_ghf == 1)||(obj.build_p == 1)
                        b = 1;
                        return;
                    elseif obj.build_p == 0
                        b = 0;
                        return;
                    end
                 end
             elseif strcmpi(l,'c') && obj.build_c == 1
                 b = 1;
                 return;
             elseif (strcmpi(l,'u') || strcmpi(l,'u0') || strcmpi(l,'a')) && obj.build_birat
                 b = 1;
                 return;
             elseif contains('wz',l,'IgnoreCase',true) && obj.build_birat
                 b = 1;
             else
                 b = 0;
             end
        end
        
        function b = check_range(obj,l,i,varargin)%l,i,j)
             %checks the range of indices for all possible variables
             %that occur in the app
             if ~isempty(varargin)
                 j = varargin{1};
             else
                 j = cell(0);
             end

             
             %check that everything has been indeed assigned
             %note that this also checks that len is appropriate
             
             if ~isempty(j)
                 
                 m = max(i,j);
                 t = min(i,j);
                 if t < 1 || m > obj.n
                     b = 0;
                     return;
                 end
                 
                 if strcmpi(l,'g') && (i >= j)
                     b = 1;
                 elseif strcmpi(l,'f') && (i+j <= obj.n-1)
                     b = 1;
                 elseif strcmpi(l,'h') && (i <= j)
                     b = 1;
                 elseif (strncmpi(l,'p',1)) && (i+j <= obj.n)
                     b = 1;
                 elseif contains('xywzu',l,'IgnoreCase',true)
                     b = 1;
                 elseif strcmpi(l,'u')||strcmpi(l,'a')||strcmpi(l,'u0')
                     b = 1;
                 else
                     b = 0;
                 end
             elseif strcmpi(l,'c') && (i > -1) && (i < obj.n+1)
                 b = 1;
                 return;
             else
                 b = 0;
             end
             
        end
        
        function b = check_mut_range(obj,type)
            %only checks that the variable of the given type
            %is on the quiver; it doesn't check existence nor indices
           if contains('gfh',type,'IgnoreCase',true)
               b = 1;
           elseif strncmpi('p',type,1)
               b = 1;
           else 
               b = 0;
           end
        end
        
        function b = check_show_range(obj,type)
           if contains('gfh',type,'IgnoreCase',true)
               b = 1;
           elseif strncmpi('p',type,1)
               b = 1;
           elseif strcmpi('u',type)
               b = 1;
           elseif strcmpi('u0',type)
               b = 1;
           elseif strcmpi('c',type)
               b = 1;
           elseif strcmpi('a',type)
               b = 1;
           else
               b = 0;
           end
        end
        
        function b = check_brack_range(obj,type)
           if contains('gfhcxy',type,'IgnoreCase',true)
               b = 1;
           elseif strncmpi('p',type,1)
               b = 1;
           else
               b = 0;
           end
        end
        
        function b = check_brack_range_WZ(obj,type)
           if contains('gfhwzuca',type,'IgnoreCase',true)
               b = 1;
           elseif strncmpi('p',type,1)
               b = 1;
           elseif strcmpi('u0',type)
               b = 1;
           else
               b = 0;
           end
        end
        
        function output = check(obj,case_type,varargin)
             len = length(varargin);
             l = cell(0); i = cell(0); j = cell(0);
             
             %attempt to pick the name and indices for a variable
             for k = 1:len
                if ischar(varargin{k})
                    l = varargin{k};
                elseif isempty(i) && isnumeric(varargin{k})
                    i = varargin{k};
                elseif isempty(j) && isnumeric(varargin{k})
                    j = varargin{k};
                end
                if ~(isempty(l) || isempty(i) || isempty(j))
                   break; 
                end
             end          
             
             if isempty(l)
                 output = {'The type of the variable has not been indicated.'};
                 return;
             elseif isempty(i)
                 output = {'The index of the variable has not been specified.'};
                 return;
             end
             
            output = cell(0);
            case_mutations = 0;
            case_show = 0;
            case_bracket = 0;
            case_bracket_WZ = 0;
            case_quiver = 0;
            %different operations have a different range of variables
            if strncmpi(case_type,'Mut',3)
                case_mutations = 1;
            elseif strncmpi(case_type,'Show',4)
                case_show = 1;
            elseif strncmpi(case_type, 'Brack',5)||strcmpi(case_type,'XY')
                case_bracket = 1;
            elseif strncmpi(case_type, 'Quiver',4)
                case_quiver = 1; 
            elseif strcmpi(case_type, 'WZ')
                case_bracket_WZ = 1;
            end
            if case_mutations || case_quiver
                if ~obj.check_mut_range(l)
                    output = {'The variable type does not belong to the quiver.'};
                    return;
                end
                if ~obj.check_mut_exception(l,i,j)
                   output = {'To mutate at φ(1,1), one needs to construct both φ- and c-functions.'};
                   return;
                end
            elseif case_show
                if ~obj.check_show_range(l)
                   output = {'The variable is of the type that can''t be displayed.'};
                   return;
                end
            elseif case_bracket
                if ~obj.check_brack_range(l)
                   output = {'The bracket can''t be taken for this type of variables.'};
                   return;
                end
            elseif case_bracket_WZ
                if ~obj.check_brack_range_WZ(l)
                   output = {'The bracket can''t be taken for this type of variables.'};
                   return;
                end
            end
            if ~obj.check_existence(l,i,j) && ~case_quiver
               output = {'The indicated variable has not been constructed.'};
               return;
            end
            if ~obj.check_range(l,i,j)
                output = {'The indicated indices are out of range.'};
                return;
            end
            
        end
        function b = check_mut_exception(obj,l,i,j)
           %if c-functions are not constructed, then we can't
           %mutate at phi(1,1)
           if isempty(l)||isempty(i)||isempty(j)
               b = 1;
               return;
           end
           if strcmpi(l,'p')&& (i == 1) && (j == 1) 
               if ~(obj.build_p && obj.build_c)
                   b = 0;
                   return;
               else
                   b = 1;
                   return;
               end
           else
               b = 1;
           end
        end
    end
    methods(Static)
        function output = BD_check(varargin)
            len = length(varargin);
            %I should have not used varargin, for I need all variables
            for i = 1:len
               if strcmpi(varargin{i},'G1_r')
                   G1_r = error_handler.extract_G(varargin{i+1});
               elseif strcmpi(varargin{i},'G2_r')
                   G2_r = error_handler.extract_G(varargin{i+1});
               elseif strcmpi(varargin{i},'G1_c')
                   G1_c = error_handler.extract_G(varargin{i+1});
               elseif strcmpi(varargin{i},'G2_c')
                   G2_c = error_handler.extract_G(varargin{i+1});
               elseif strcmpi(varargin{i},'n')
                   n = varargin{i+1};
               elseif strcmpi(varargin{i},'build_birat')
                   build_birat = varargin{i+1}; 
               end
            end
            if n < 3
               output = cell(1);
               output{1} = 'The value of n must be bigger than or equal to 3.';
               return;
            end
            b = error_handler.BD_cor(G1_r,G2_r,n);
            if b
                b = error_handler.BD_cor(G1_c,G2_c,n);
            else
                output = cell(1);
                output{1} = 'An error in the row BD-data.';
                return;
            end
            if ~b
               output = cell(1);
               output{1} = 'An error in the column BD-data.';
               return;
            end
            b = error_handler.BD_acyclicity(G1_r,G2_r,G1_c,G2_c,n);
            if ~b
               output = cell(1);
               output{1} = 'The BD-graph contains a cycle.';
               output{1} = append(output{1}, ' The cluster structure can''t be constructed.');
               return; 
            end
            
            if build_birat && G1_r(1) == 0
                output = cell(1);
                output{1} = 'To construct U, need a non-empty row BD-data.';
                return;
            end
            
            output = cell(0);
        end
        function b = BD_cor(G1,G2,n)
            %returns 1 if the data is good
            len_G1 = length(G1);
            len_G2 = length(G2);
            %here we basically check that the BD map is a bijection on the
            %right domain and range
            if len_G1 ~= len_G2
                b = 0;
                return;
            end
            if (len_G1 ~= length(unique(G1)))|| (len_G2 ~= length(unique(G2)))
                b = 0;
                return;
            end
            if len_G1 >n-1 || len_G2 > n-1
                b = 0;
                return;
            end
            
            %now a more subtle one: check the map is orthogonal and
            %oriented
            bd = BD_data(G1,G2,n);
            runs1 = bd.runs(G1);

            %the same number of runs we have at this point!
            %check the runs correspond to the runs with the same length
            for i = 1: length(runs1(:,1))
                rlen = runs1(i,2) - runs1(i,1);
                if rlen == 0
                    continue;
                end
                [x,y] = bd.ipm(bd.G2_r,bd.T_r(runs1(i,1)));
                x = x+1; %to get the actual run;
                if y-x ~= rlen
                    b = 0;
                    return;
                end
                
                %along the lines, check that BD-data is oriented
                for j = runs1(i,1) : runs1(i,2)-2
                   if bd.T_r(j+1) ~= bd.T_r(j) + 1
                      b = 0;
                      return;
                   end
                end
            end
            
            %finally, if nothing has been caught, set b = 1:
            b = 1;
        end
        function b = BD_acyclicity(G1_r,G2_r,G1_c,G2_c,n)
           if nargin == 3
               n = G1_c;
               G1_c = G1_r;
               G2_c = G2_r;
           end
           bd = BD_data(G1_r,G2_r,G1_c,G2_c,n);
           cc = bd.max_paths();
           if isempty(cc)
               b = 0;
               return;
           end
           %now check for the cycles that we couldn't catch
           %the sum is a nice invariant
           [N,~] = size(cc);
            for i = 1:N
                if cc(i,1) == -1
                    cc(i,1) = 0;
                end
                if cc(i,2) == -1
                    cc(i,2) = 0;
                end
            end
            s = sum(sum(cc));
            if s ~= 2*n*(n-1)
                b = 0;
            else
                b = 1;
            end
           %I'll use a standard Matlab's function ...
           %Also, I assume that the BD_data has passed BD_cor test already
           %create digraph by specifying the edges
           
           %here's my algorithm: leave a coin at each node
           %then attempt to pick up all coins but not in a cycle...
            %...

        end
       function res = extract_G(str)
           %converts a string to a numeric BD-array
           G = extract(str,digitsPattern);
           if isempty(G)
               res = 0;
               return;
           else
               res = cellfun(@str2num, G);
           end
       end
       
    end
end