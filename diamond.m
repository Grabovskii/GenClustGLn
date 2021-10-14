classdef diamond < handle 
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    
    %just to note somewhere: I verified on paper that 
    %f_{n-l,l} = phi_{n-l,l}; i.e., no negative signs involved
    properties
        n
        F, P %cell arrays of matrices of which we then take determinants
        
        build_p_
    end
    properties(Constant)
       %default settings
       BUILD_P = 0; %for large quivers it's impossible to construct them in a timely manner 
       OPTIMIZATION = 1; %the lower row of phi's from the quiver is defined in
       %the equivalent way as phi_{n-l,l} := f_{n-l,l}. helps
       %computationally
    end
    
    methods
        function obj = diamond(n,varargin)
            X = sym('X', [n n]);
            Y = sym('Y', [n n]);
            obj.n = n;
            obj.build_p_ = obj.BUILD_P; %launch first default
            obj.F = cell(n);
            for k=1:n-1%n-2
                for l = 1:n-k%n-k-1
                    A = [X(:,n-k+1:n) Y(:,n-l+1:n)];
                    obj.F{k,l} = A(n-k-l+1:n,:);
                end
            end
            
            
            %A very slow way to define phi's
            %Probably for higher n's I need to find a different way
            %to define phi's....
            obj.P = cell(n);
            if obj.OPTIMIZATION
                for k =1:obj.n-1
                    obj.P{k,n-k} = obj.F{k,n-k};
                end
            end
            
            if length(varargin) > 1
               if strcmpi(varargin{1},'build_p')
                  obj.build_p_ = varargin{2}; 
               end
            end
            
            if obj.build_p_ > 0
                U = inv(X)*Y; %takes time ...
                Ucells = cell(n-1); %powers of U >= 1 till n-k-l+1
                Ucells{1} = U;
                %Ucells{1} = zeros(n,n-k+1,n); %U^0
                for i = 2:(n-1)
                   Ucells{i} = Ucells{i-1} * U; %takes a loooot of time
                end
                I = poisson_bracket.sym_eye(n);
                %initialize phi's (U^0 and U^1);
                for k = 1:(n-1)
                    for l = 1:(n-k)
                        if (obj.OPTIMIZATION)&&(k+l == n)
                           continue;
                        end
                        obj.P{k,l} = I;
                        obj.P{k,l}(:,1:k) = I(:,n-k+1:n);
                        obj.P{k,l}(:,k+1:k+l) = U(:,n-l+1:n);
                    end
                end
                %now correct them with U's
                for k=1:(n-1)
                    for l = 1:(n-k)
                        if (obj.OPTIMIZATION)&&(k+l == n)
                           continue;
                        end
                        for i =2:(n-k-l+1)
                            obj.P{k,l}(:,k+l-1+i:k+l-1+i) = Ucells{i}(:,n:n);
                        end
                    end
                end
            end
        end
        function show_f(obj)
            for k=1:obj.n-1%n-2
                for l = 1:obj.n-k%n-k-1
                    fprintf('F(%d,%d) = \n',k,l);
                    diamond.fpm(obj.F{k,l});
                    fprintf('\n');
                end
            end
            
        end
        function h = f(obj,i,j)
            h = det(obj.F{i,j});
        end
        function p = phi(obj,k,l,c)
            %looks cumbersome, I agree
            %the parameter c determines whether we multiply by X or W
            %where the W comes from tests with the birational map
            if nargin == 3
                c = 0;
            end
            p = det(obj.P{k,l});
            if obj.OPTIMIZATION && (k+l == obj.n)
                return;
            end
            if (nargin == 3)
                X = sym('X', [obj.n obj.n]);
                p = det(X)^(obj.n-k-l+1) * p;
            elseif c == 1
                W = sym('W', [obj.n obj.n]);
                p = det(W)^(obj.n-k-l+1) * p;
            end
            if mod(obj.n,2) == 0
               p = p*(-1)^(k*(l+1));
            else
               p = p* (-1)^((obj.n-1)/2 + k*(k-1)/2 + l*(l-1)/2);
            end
            p = simplifyFraction(p,'Expand',true);
        end
    end
    methods(Static)
        function fpm(Q)
            %prints the matrix Q
            [m,~] = size(Q);
            for i = 1:m
                fprintf('|');
               for j = 1:m
                  fprintf('%s',char(Q(i,j))); 
                  if j < m
                     fprintf(' '); 
                  else
                     fprintf('|\n'); 
                  end
               end
            end
        end 
    end
end

