classdef experiment < double_full
    properties
        G1_r,G2_r,G1_c,G2_c,n %BD pair
    end
    
    methods
        function obj = experiment(G1_r,G2_r,G1_c,G2_c,n)  
            if nargin == 3
                n = G1_c;
                G1_c = G1_r;
                G2_c = G2_r;
            end
            obj@double_full(G1_r,G2_r,G1_c,G2_c,n);
            obj.G1_c = G1_c;
            obj.G2_c = G2_c;
            obj.G1_r = G1_r;
            obj.G2_r = G2_r;
            obj.n = n;
        end
        
%         function result = borelPlus(obj,type,i,j)
%            %returns conditions that will make the left derivative belong to the upper-triangular matrices 
%             f = obj.get_f(type,i,j);
%             lx = obj.b.gradXL(f);
%             ly = obj.b.gradYL(f);
%             sx = 0;
%             sy = 0;
%             for i = 1:obj.n^2
%                if lx(i) ~= sym('0')
%                    sx = sx + 1;
%                end
%                if ly(i) ~= sym('0')
%                   sy = sy+1; 
%                end
%             end
%             v = sym(zeros(1,sx+sy)); %vector of nonzero derivatives of f
%             d = sym(zeros(1,sx+sy)); %here we store what variables we've taken derivatives of
%             s = 1;
%         end

      function result = infsSymmetry(obj,type,i,j)
            choice = 'LL';
            %choice = 'LL';
            %returns the infinitesimal relations for semi-invariance
            %i.e., invariance that comes with a character
            f = obj.get_f(type,i,j);
            if strcmp(choice,'LL')
                gx = simplifyFraction(obj.b.gradXL(f),'Expand',true);
                gy = simplifyFraction(obj.b.gradYL(f),'Expand',true);
            elseif strcmp(choice,'RL')
                gx = simplifyFraction(obj.b.gradXR(f),'Expand',true);
                gy = simplifyFraction(obj.b.gradYL(f),'Expand',true);
            elseif strcmp(choice,'LR')
                gx = simplifyFraction(obj.b.gradXL(f),'Expand',true);
                gy = simplifyFraction(obj.b.gradYR(f),'Expand',true);          
            elseif strcmp(choice,'RR')
                gx = simplifyFraction(obj.b.gradXR(f),'Expand',true);
                gy = simplifyFraction(obj.b.gradYR(f),'Expand',true);
            end
            U = sym('U',obj.n);
            V = sym('V',obj.n);
            XX = sym('XX',obj.n);
            YY = sym('YY',obj.n);
            
            mf = subs(f,[obj.s.X,obj.s.Y],[XX,YY]);
            term = trace(gx*U)+trace(gy*V);
            termt = subs(term,[obj.s.X,obj.s.Y],[XX,YY]);
            row = mf*term - f*termt; %a row in the equations
            rowCoeff = sym(zeros([1,2*obj.n^2])); %coefficients at U and V
            eii = sym(zeros(obj.n));
            z = eii;
            for i = 1:obj.n^2
               eii(i) = sym('1');
               rowCoeff(i) = subs(row,[U,V],[eii,z]);
               rowCoeff(i+obj.n^2) = subs(row,[U,V],[z,eii]);
               eii(i) = sym('0');
            end
            %хотя бы 2*obj.n^2 уравнений нужно...  наверное
            num_rows = 5*obj.n^2;
            A = zeros(num_rows,2*obj.n^2);
            for i = 1:num_rows
               rx = randi([1,1000], obj.n);
               rxx = randi([1,1000],obj.n);
               ry = randi([1,1000], obj.n);
               ryy = randi([1,1000],obj.n);
               A(i,:) = subs(rowCoeff,[obj.s.X,obj.s.Y,XX,YY],[rx,ry,rxx,ryy]);
            end
            N = null(A);
            N = transpose(rref(transpose(N)));
            %clear the array
            [rN,cN] = size(N);
            for i = 1:rN*cN
               if(N(i)<1e-12)
                   N(i) = sym('0');
               end
            end
            result = sym(zeros([1,length(N(1,:))]));
            %now comes interpretation of N
            for j = 1:length(N(1,:))
                for i=1:obj.n^2
                    if N(i,j) ~= 0
                        result(j) = result(j) + U(i)*N(i,j);
                    end
                    if N(i+obj.n^2,j) ~= 0
                        result(j) = result(j) + V(i)*N(i+obj.n^2,j); 
                    end
                end
            end
            
        end

%       function result = infsSymmetry(obj,type,i,j)
%             %choice = 'RL';
%             choice = 'LL';
%             %returns the infinitesimal relations for semi-invariance
%             %i.e., invariance that comes with a character
%             f = obj.get_f(type,i,j);
%             if strcmp(choice,'LL')
%                 gx = simplifyFraction(obj.b.gradXL(f),'Expand',true);
%                 gy = simplifyFraction(obj.b.gradYL(f),'Expand',true);
%             elseif strcmp(choice,'RL')
%                 gx = simplifyFraction(obj.b.gradXR(f),'Expand',true);
%                 gy = simplifyFraction(obj.b.gradYL(f),'Expand',true);
%             elseif strcmp(choice,'LR')
%                 gx = simplifyFraction(obj.b.gradXL(f),'Expand',true);
%                 gy = simplifyFraction(obj.b.gradYR(f),'Expand',true);          
%             elseif strcmp(choice,'RR')
%                 gx = simplifyFraction(obj.b.gradXR(f),'Expand',true);
%                 gy = simplifyFraction(obj.b.gradYR(f),'Expand',true);
%             end
%             U = sym('U',obj.n);
%             V = sym('V',obj.n);
%             term = trace(gx*U)+trace(gy*V);
%             z = sym(zeros(obj.n));
%             eii = z;
%             expr = sym(zeros(obj.n,2*obj.n));
%             for i = 1:obj.n^2
%                expr(i) = -diff(f,obj.s.X(i))*term+f*diff(term,obj.s.X(i));
%                expr(i) = subs(expr(i),[obj.s.X,obj.s.Y],[z,z]);
%                expr(i+obj.n^2) = -diff(f,obj.s.Y(i))*term+f*diff(term,obj.s.Y(i));
%                expr(i+obj.n^2) = subs(expr(i+obj.n^2),[obj.s.X,obj.s.Y],[z,z]);
%             end
%             A = zeros(2*obj.n^2,2*obj.n^2);
%             for i =1:2*obj.n^2 
%                 if expr(i) == 0
%                     continue;
%                 end
%                 eii = z;
%                 for j=1:obj.n^2
%                    eii(j) = 1;
%                    A(i,j) = subs(expr(i),[U,V],[eii,z]);
%                    A(i,j+obj.n^2) = subs(expr(i),[U,V],[z,eii]);
%                    eii(j) = 0;
%                 end
%             end
%             N = null(A);
%             result = sym(zeros([1,length(N(1,:))]));
%             %now comes interpretation of N
%             for j = 1:length(N(1,:))
%                 for i=1:obj.n^2
%                     if N(i,j) ~= 0
%                         result(j) = result(j) + U(i)*N(i,j);
%                     end
%                     if N(i+obj.n^2,j) ~= 0
%                         result(j) = result(j) + V(i)*N(i+obj.n^2,j); 
%                     end
%                 end
%             end
%         end

        function result = infSymmetry(obj,type,i,j)
            %returns infinitesimal relations that yield the derivatives
            %equal to zero
            f = obj.get_f(type,i,j);
            lx = simplifyFraction(obj.b.gradXL(f),'Expand',true);
            ly = simplifyFraction(obj.b.gradYL(f),'Expand',true);
            %how many non-zero derivatives do we have?
            sx = 0;
            sy = 0;
            for i = 1:obj.n^2
               if lx(i) ~= sym('0')
                   sx = sx + 1;
               end
               if ly(i) ~= sym('0')
                  sy = sy+1; 
               end
            end
            v = sym(zeros(1,sx+sy)); %vector of nonzero derivatives of f
            d = sym(zeros(1,sx+sy)); %here we store what variables we've taken derivatives of
            s = 1;
            for i = 1:obj.n^2
                if lx(i) ~= sym('0')
                    v(s) = lx(i);
                    d(s) = obj.s.X(i);
                    s = s +1;
                end
                if ly(i) ~= sym('0')
                    v(s) = ly(i);
                    d(s) = obj.s.Y(i);
                    s = s +1 ;
                end
            end
            A = sym(zeros(sx+sy));
            for i = 1:sx+sy
                rx = randi([1,1000], obj.n);
                ry = randi([1,1000], obj.n);
                A(i,:) = subs(v,[obj.s.X,obj.s.Y],[rx,ry]);
            end
            ker = null(A);
            dim = length(ker(1,:));
            result = sym(zeros(1,dim));
            for j = 1:dim
               for i = 1:sx+sy
                  if ker(i,j)~=0
                      result(j) = result(j) + ker(i,j)*d(i);
                  end
               end
            end
        end
        
       function B = rpattern(obj,A)
           B = subs(A,[obj.s.X obj.s.Y], [randi([1,100], [obj.n obj.n]) randi([1,100],[obj.n obj.n])]);
           for i = 1:obj.n
               for j = 1:obj.n
                if B(i,j) ~= 0
                    B(i,j) = 1;
                end
               end
           end
        end
        
        function B = spattern(obj,A)
           %a more strict check
           %очень полезная функция!!!
           B1 = obj.rpattern(A);
           B2 = obj.rpattern(A);
           if B1== B2
               B = B1;
               return;
           else 
              B = obj.spattern(A);
           end
        end
 
    end
end

