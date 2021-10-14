classdef exphi < BD_data

    
    properties 
       s %the shores
       P %the Phi's without pluggging in U = X^{-1} Y
       b %Poisson bracket
       U
       U_s %here I compute only once U_spec = X^{-1} Y and store
    end
    
    methods
        function obj = exphi(G1_r,G2_r,G1_c,G2_c,n)
            if nargin == 3
               n = G1_c;
               G1_c = G1_r;
               G2_c = G2_r;
            end
           obj@BD_data(G1_r,G2_r,G1_c,G2_c,n);
           obj.s = shores(G1_r,G2_r,G1_c,G2_c,n); 
           obj.b = poisson_bracket(G1_r,G2_r,G1_c,G2_c,n);
           
           %create Phi's
            obj.U = sym('U', [obj.n, obj.n]);
            U = obj.U;
            obj.U_s = subs(U,U,inv(obj.b.X)*obj.b.Y);
            Ucells = cell(n-1); %powers of U >= 1 till n-k-l+1
            Ucells{1} = U;
            %Ucells{1} = zeros(n,n-k+1,n); %U^0
            for i = 2:(n-1)
               Ucells{i} = Ucells{i-1} * U; %takes a loooot of time
            end
            obj.P = cell(n);
            I = sym('I',[n n]); %symbolic identity matrix
            for i = 1:n
                for j = 1:n
                    if i~=j
                        I(i,j) = sym('0');
                    else
                        I(i,j) = sym('1');
                    end
                end
            end
            %initialize phi's (U^0 and U^1);
            for k = 1:(n-1)
                for l = 1:(n-k)
                    obj.P{k,l} = I;
                    obj.P{k,l}(:,1:k) = I(:,n-k+1:n);
                    obj.P{k,l}(:,k+1:k+l) = U(:,n-l+1:n);
                end
            end
            %now correct them with U's
            for k=1:(n-1)
                for l = 1:(n-k)
                    for i =2:(n-k-l+1)
                        obj.P{k,l}(:,k+l-1+i:k+l-1+i) = Ucells{i}(:,n:n);
                    end
                end
            end
        end
        function result = lterm(obj,p,g)
%             com = poisson_bracket.comm(obj.gradU(p),obj.U);
%             gr = obj.b.EL(g);%obj.b.gradX(g)*obj.b.X;
%             result = poisson_bracket.sym_trace(obj.b.R_plus45(com)*gr);
            com = poisson_bracket.comm(obj.gradU(p),obj.U);
            gr = obj.b.EL(g);%obj.b.gradX(g)*obj.b.X;
            result = poisson_bracket.sym_trace(obj.b.Tinv(obj.b.TTinv((com)))*gr);
        end
        function result = ltest(obj,i,j,k,l)
           result = obj.lterm(log(det(obj.P{i,j})),log(det(obj.s.G{k,l}))); 
           result = subs(result,obj.U,inv(obj.b.X)*obj.b.Y);
           result = simplifyFraction(result,'Expand',true);
        end
        function D = gradU(obj,p)
            D = poisson_bracket.sym_zeros(obj.n);
            for i =1:obj.n
                for j =1:obj.n
                    D(i,j) = diff(p,obj.U(j,i));
                end
            end
        end
        function display_grads(obj)
            k = 2; l = 1;
           for i = 1:obj.n-1
               for j =1:obj.n-i
                   %for the very gradients:
                   %fprintf('gradU(detP{%d,%d}) = \n',i,j);
                   %simplifyFraction(obj.gradU(det(obj.P{i,j})),'Expand',true)
                   %for the commutators of gradients with U:
                   %fprintf('[gradU(detP{%d,%d}),U] = \n',i,j);
%                    if (i~= j || i ~= 1)
%                        fprintf('R([gradU(detP{%d,%d}),U]) = \n',i,j);
%                        result = poisson_bracket.comm(obj.gradU(det(obj.P{i,j})),obj.U);
%                        result = simplifyFraction(obj.b.R_plus45(result),'Expand',true);
%                        r = poisson_bracket.pattern(result)
%                    end
                    if (i~=j || i~=1)
                       fprintf('(T*/(1-T*)([gradU(log(P(%d,%d)),U]),log(g%d%d)) = \n',i,j,k,l);
                       result = obj.ltest(i,j,k,l)
                    end
               end
           end
        end
%         function r = eqq(obj,k,l,i,j)
%            Ekl = poisson_bracket.e(k,l,obj.n);
%            r = poisson_bracket.sym_trace(poisson_bracket.comm(obj.U,obj.b.R_plus45star(Ekl))*obj.gradU(det(obj.P{i,j})));
%            r = simplifyFraction(r,'Expand',true);
%         end
    end
end