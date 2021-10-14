classdef exotic
    %this class produces a family of functions from Misha's
    %article on Exotic cluster structures on SLn
    %upd: I guess I wrote this class only to check that it coincides
    %with special cases from plethora; it indeed does

    
    properties (Access = public)
        Theta, Phi, Psi %the families of functions
        n
        X,Y
    end
    properties (Access = private)
        U
    end
    
    methods (Access = public)
        function obj = exotic(n)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.n = n;
            obj.X = sym('X',[n n]);
            obj.Y = sym('Y',[n n]);
            X_curly = exotic.sym_zeros(n-1,n+1);
            Y_curly = exotic.sym_zeros(n-1,n+1);
            X_curly(:,1:n) = obj.X(2:n,:);
            Y_curly(:,2:n+1) = obj.Y(1:n-1,:); 
            k = floor((n+1)/2); 
            obj.U = exotic.sym_zeros(k*(n-1),(k+1)*(n+1));
            for i=1:k
               obj.U((i-1)*(n-1)+1:i*(n-1),(i-1)*(n+1)+1:i*(n+1)) = Y_curly;
               obj.U((i-1)*(n-1)+1:i*(n-1),i*(n+1)+1:(i+1)*(n+1)) = X_curly;
            end
            obj.Theta = obj.ini_theta();
            obj.Psi = obj.ini_psi();
            obj.Phi = obj.ini_phi();
        end
        function p = phi(obj,i)
           p = det(obj.Phi{i});
        end
        function p = psi(obj,i)
            p = det(obj.Psi{i});
        end
        function p = theta(obj,i)
            p = det(obj.Theta{i});
        end
        function display_fns(obj)
           n1 = size(obj.Theta);
           n2 = size(obj.Psi);
           n3 = size(obj.Phi);
           for i =1:n1
               fprintf('Theta(%d) = ',i);
               obj.Theta{i}
           end
           for i = 1:n2
              fprintf('Psi(%d) = ',i);
              obj.Psi{i}
           end
           for i = 1:n3
              fprintf('Phi(%d) = ',i);
              obj.Phi{i}
           end
        end
    end
    methods(Access = private)
        function f = ini_theta(obj)
            f = cell(obj.n-1);
            for i = 1:obj.n-1
                f{i} = obj.X(obj.n-i+1:obj.n,obj.n-i+1:obj.n);
            end
        end
        
        function f = ini_phi(obj)
            k = floor((obj.n+1)/2);
            N = k*(obj.n-1);
            f = cell(N);
            for i = 1:N
               f{i} = obj.U(N-i+1:N,k*(obj.n+1)-i+1:k*(obj.n+1));
            end
        end
        
        function f = ini_psi(obj)
            k = floor((obj.n+1)/2);
            N = k*(obj.n-1);
           if mod(obj.n,2) == 0
              M = N; 
           else
              M = N-obj.n+1; 
           end
           f = cell(M);
           for i =1:M
              f{i} = obj.U(N-i+1:N, k*(obj.n+1)-i+2:k*(obj.n+1)+1);
           end
        end
    end
    methods (Static)
        function Z = sym_zeros(m,n)
           Z = sym('Z', [m n]);
           for i = 1:m 
               for j = 1:n
                    Z(i,j) = sym('0');
               end
           end
        end
    end
end

