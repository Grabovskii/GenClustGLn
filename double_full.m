classdef double_full < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
       TECH = 0;%1; %if true, we get additional methods of manipulating individual blocks; 
                 %need for testing Misha's formulae only.
    end
    
    properties (Access = public)
        d, s, c %functions: diamond, shores, and isolated c's
        q %the quiver for the double
        b %the bracket
        build_ghf_
        build_c_
        build_p_
        bracket_
    end
   
    properties (Constant)
        %Default settings
       BUILD_GHF = 1; %construct or not g-,h-,f- functions
       BUILD_P = 0; %construct the strict upper triangle with phi or not
       BUILD_C = 0; %the isolated vertices
       BRACKET = 1; %construct the Poisson bracket as well
    end
    
    methods
        function obj = double_full(G1,G2,varargin)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            lvar = length(varargin);
            %default values
            obj.build_ghf_ = obj.BUILD_GHF;
            obj.build_c_ = obj.BUILD_C;
            obj.build_p_ = obj.BUILD_P;
            obj.bracket_ = obj.BRACKET;
            
            if lvar > 0
                if lvar == 1 
                    n = varargin{1};
                    G1_r = G1;
                    G2_r = G2;
                    G1_c = G1_r;
                    G2_c = G2_r;
                elseif lvar > 1 && ischar(varargin{2})
                    n = varargin{1};
                    G1_r = G1;
                    G2_r = G2;
                    G1_c = G1_r;
                    G2_c = G2_r;     
                end
            end
            if lvar > 2
               if ~(ischar(varargin{1}) || ischar(varargin{2}) || ischar(varargin{3}))
                   G1_r = G1;
                   G2_r = G2;
                   G1_c = varargin{1};
                   G2_c = varargin{2};
                   n = varargin{3};
               end
            end
            
            if lvar > 1
                for i = 1:lvar
                    if ischar(varargin{i})
                        str = varargin{i};
                        if strcmpi(str,'build_ghf')
                            obj.build_ghf_ = varargin{i+1};
                        elseif strcmpi(str,'build_c')
                            obj.build_c_ = varargin{i+1};
                        elseif strcmpi(str,'bracket')
                            obj.bracket_ = varargin{i+1};
                        elseif strcmpi(str, 'build_p')||strcmpi(str,'build_phi')
                            obj.build_p_ = varargin{i+1};
                        end
                    end
                end
            end

            if obj.build_ghf_
                if ~obj.TECH
                    obj.s = shores(G1_r,G2_r,G1_c,G2_c,n);
                elseif obj.TECH
                    obj.s = shorestech(G1_r,G2_r,G1_c,G2_c,n);
                end
                obj.d = diamond(n,'build_p',obj.build_p_);
            end
            if obj.build_c_ == 1
               obj.c = cisol(n); 
            end
            if obj.bracket_ == 1
               obj.b = poisson_bracket(G1_r,G2_r,G1_c,G2_c,n); 
            end
            obj.q = double_gln(G1_r,G2_r,G1_c,G2_c,n);
        end
        
        function y = get_y(obj,type,i,j)
           %extract the y-coordinate of type(i,j) (no logarithm is taken)
           y = 1;
           ind = obj.q.ind(i,j,type);
           for i = 1:obj.q.dim
                if obj.q.adj_matrix(ind,i) ~= 0
                    y = y*(obj.get_f(i)^obj.q.adj_matrix(ind,i));
                end
           end
        end
        
        
        function f = get_f(obj,type,i,j)
           %this function returns a function of type
           %that might be 'g', 'phi', 'h', or 'f'
           %of index (i,j)
           if nargin == 4
               if strcmp(type,'g')
                  f = obj.s.g(i,j); 
               elseif strcmp(type,'f')
                  f = obj.d.f(i,j);
               elseif strcmp(type,'h')
                  f = obj.s.h(i,j);
               elseif strcmp(type,'phi')||strcmp(type,'p')
                   f = obj.d.phi(i,j);
               end
           elseif nargin == 2
               iind = obj.q.back_d(type);%wierd, but type now is index in adj_matrix 
               f = obj.get_f(iind{1},iind{2},iind{3});
           end
        end
        
        function f = get_m(obj,type,i,j)
           %returns a matrix that corresponds
           %to the function
           if nargin == 4
               if strcmp(type,'g')
                  f = obj.s.G{i,j}; 
               elseif strcmp(type,'f')
                  f = obj.d.F{i,j};
               elseif strcmp(type,'h')
                  f = obj.s.H{i,j};
               elseif strcmp(type,'phi')||strcmp(type,'p')
                   f = obj.d.P{i,j};
               end
           elseif nargin == 2
               iind = obj.q.back_d(type);%wierd, but type now is index in adj_matrix 
               f = obj.get_f(iind{1},iind{2},iind{3});
           end
        end
        
        function mutate_B(obj,s)
            n = obj.q.n;
            for i=n:-1:s
               fprintf('h"(%d,%d) = \n',s,i);
               disp(obj.mutate('h',s,i));
            end
            for i = 1:n-s
               fprintf('f"(%d,%d) = \n',i,n-s-i+1);
               disp(obj.mutate('f',i,n-s-i+1));
            end
            for i = s:-1:2
                fprintf('g"(%d,%d) = \n',s,i);
                disp(obj.mutate('g',s,i));
            end
        end
        
        function mutate_Bi(obj,s)
            n = obj.q.n;
            %just a test
            %that mutate_B in the opposite direction
            for i = 2:s%s:-1:2
                fprintf('g"(%d,%d) = \n',s,i);
                disp(obj.mutate('g',s,i));
                %fprintf('\n');
            end
            for i = n-s:-1:1
               fprintf('f"(%d,%d) = \n',i,n-s-i+1);
               disp(obj.mutate('f',i,n-s-i+1));
               %fprintf('\n');
            end
            for i=s:n
                fprintf('h"(%d,%d) = \n', s,i);
               disp(obj.mutate('h',s,i));
               %fprintf('\n');
            end
        end
        
        function result = mutate(obj,type,i,j)

            %this function mutates at type(i,j), where type
            %can be 'phi', 'g', 'f', 'h'
            %returns the mutated function and does not affect the
            %properties of the object
            if obj.q.isfrozen(obj.q.ind(i,j,type))
                result = NaN;
                fprintf('This vertex is frozen!');
                return;
            end
            except = ((i == 1) && (j == 1) && (strcmp(type,'phi')||strcmp(type,'p')));
            N = obj.q.dim;
            l = sym('1'); r = sym('1'); %left and right monomial in the mutation
            ind = obj.q.ind(i,j,type);
            if ind == 0
               fprintf('There is no function to mutate.\n');
               result = sym('0');
               return; 
            end
            if ~except
                for k = 1:N
                   if obj.q.adj_matrix(ind,k) > 0
                      l = l * obj.get_f(k);
                   elseif obj.q.adj_matrix(ind,k) < 0
                      r = r * obj.get_f(k);
                   end
                end
                result = (l + r)/(obj.get_f(ind));
                result = simplifyFraction(result,'Expand',true);           
            else
                u = obj.get_f('phi',2,1);
                v = obj.get_f('phi',1,2);
                result = obj.c.get_c(0)*(v^obj.q.n);
                for r = 1:obj.q.n
                   result = result + obj.c.get_c(r) * (u^r) * v^(obj.q.n-r);
                end
                result = result/obj.get_f('phi',1,1);
                result = simplifyFraction(result,'Expand',true);
            end
            
            % ALTER FUNCTIONS AND THE ADJOINT MATRIX IF NECESSARY
            if strcmp(type,'g')
                %obj.s.set_g(result,i,j);%
                obj.s.G{i,j} = result;
            elseif strcmp(type,'h')
                %obj.s.set_h(result,i,j);
                obj.s.H{i,j} = result;
            elseif strcmp(type,'phi') || strcmp(type,'p')
                obj.d.P{i,j} = result;
            elseif strcmp(type,'f')
                obj.d.F{i,j} = result;
            end
            %update the adjoint matrix
            obj.q.mutate(type,i,j);
        end
        function show_nbd(obj,type,i,j)
           say = strcat('The neighborhoud of ',type',string(i),string(j),'\n');
           fprintf(say);
           index =  obj.q.ind(i,j,type);
           for k = 1:obj.q.dim
               if obj.q.adj_matrix(index,k) ~= 0
                   t = obj.q.back_d(k);
                   ss = strcat(t{1},string(t{2}),'_',string(t{3}),':\n');
                   fprintf(ss);
                   obj.get_m(t{1},t{2},t{3})
               end
           end
        end
        function subs_vars(obj,W,Z,diamond_val)
           %this function substitutes everywhere X with W
           %and Y with Z. Needed for testbirat (computing the birational
           %map between two cluster structures)
           if nargin == 1
               W = sym('W',obj.s.n);
               Z = sym('Z',obj.s.n);
           end
           if obj.BRACKET
              obj.b.XY_to_WZ; 
           end
           for i = 1:obj.s.n
               for j = 1:i
                    obj.s.G{i,j} = subs(obj.s.G{i,j},obj.s.X,W);
                    obj.s.G{i,j} = subs(obj.s.G{i,j},obj.s.Y,Z);
               end
           end
           for j = 1:obj.s.n
              for i = 1:j 
                  obj.s.H{i,j} = subs(obj.s.H{i,j},obj.s.X,W);
                  obj.s.H{i,j} = subs(obj.s.H{i,j},obj.s.Y,Z);
              end
           end
           if nargin==4 
               if diamond_val == 1
                   for i =1:(obj.s.n-2)
                      for j = 1:(obj.s.n-i-1)
                         obj.d.F{i,j} = subs(obj.d.F{i,j},obj.s.X,W); 
                         obj.d.F{i,j} = subs(obj.d.F{i,j},obj.s.Y,Z);
                      end
                   end
                   for i =1:(obj.s.n-1)
                      for j = 1:(obj.s.n-i)
                         obj.d.P{i,j} = subs(obj.d.P{i,j},obj.s.X,W); 
                         obj.d.P{i,j} = subs(obj.d.P{i,j},obj.s.Y,Z);
                      end
                   end
                   for i=1:(obj.s.n+1)
                      obj.c.C{i} = subs(obj.c.C{i},obj.s.X,W); 
                      obj.c.C{i} = subs(obj.c.C{i},obj.s.Y,Z);
                   end
               end
           end
        end
        
        function subs_diag(obj,Q)
           W = sym('W',obj.s.n);
           for i = 1:obj.s.n
               for j = 1:i
                    obj.s.G{i,j} = subs(obj.s.G{i,j},obj.s.X,Q);
                    obj.s.G{i,j} = subs(obj.s.G{i,j},obj.s.Y,Q);
                    obj.s.G{i,j} = subs(obj.s.G{i,j},W,Q);
               end
           end
           for j = 1:obj.s.n
              for i = 1:j 
                  obj.s.H{i,j} = subs(obj.s.H{i,j},obj.s.X,Q);
                  obj.s.H{i,j} = subs(obj.s.H{i,j},obj.s.Y,Q);
                  obj.s.H{i,j} = subs(obj.s.H{i,j},W,Q);
              end
           end
        end
    end
end

