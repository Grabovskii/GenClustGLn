classdef extended_T < BD_data
    %this class realizes methods that allow to work with T (or gamma in
    %Misha's notation) on the level of Lie algebras
    %I pulled these methods from poisson_bracket.m to make the latter file
    %easier to read and more compact
    
    properties
        G1_pattern_r, G2_pattern_r %the way the projection from gln to glG1 works
        G1_pattern_c, G2_pattern_c %the way the projection from gln to glG1 works
    end
    
    methods
        function obj = extended_T(G1_r,G2_r,G1_c,G2_c,n)
            if nargin == 3
                n = G1_c;
                G1_c = G1_r;
                G2_c = G2_r;
            end
            obj@BD_data(G1_r,G2_r,G1_c,G2_c,n);
            obj.G1_pattern_r = obj.ini_G_pattern(G1_r,n);
            obj.G2_pattern_r = obj.ini_G_pattern(G2_r,n);
            obj.G1_pattern_c = obj.ini_G_pattern(G1_c,n);
            obj.G2_pattern_c = obj.ini_G_pattern(G2_c,n);
        end
        
        function h = cartan_coeff(obj,x)
            type = strcmp(class(x),'sym');
            if type == 0
                h = zeros(1,obj.n-1); %coefficients of the Cartan part in terms of simple root vectors
                h(1) = x(1,1);
                for i=2:obj.n-1
                    h(i) = x(i,i) + h(i-1);
                end
            else 
                h = sym('h', [obj.n-1 1]);
                h(1,1) = x(1,1);
                for i=2:obj.n-1
                    h(i) = x(i,i) + h(i-1);
                end
            end 
        end
        
        function y = expTr(obj,x)
            %this is not the real exponent and desgined solely
            %to act properly upon upper- and lower- triangular matrices
            %outside the diagonal, it acts as Tr; 
            %the diagonal is fixed
           type = strcmp(class(x),'sym');
           if type
              d = poisson_bracket.sym_eye(obj.n).*x;
           else
              d = eye(obj.n).*x;
           end 
           y = obj.Tr(x-d)+d;
        end
        
        function y = Tr(obj,x)
           y = obj.T(x,0); 
        end
        function y = Tc(obj,x)
           y = obj.T(x,1); 
        end
        function y = Tinvr(obj,x)
           y = obj.Tinv(x,0); 
        end    
        function y = Tinvc(obj,x)
           y = obj.Tinv(x,1); 
        end
        
        function y = T(obj,x,rc)
            %if rc == 0, this is for gamma^r
            %if rc == 1, works for gamma^c
            type = strcmp(class(x),'sym');
            x = poisson_bracket.proj_sln(x);
            h = obj.cartan_coeff(x);
            if rc == 0
               G1_pattern = obj.G1_pattern_r;
            elseif rc == 1
               G1_pattern = obj.G1_pattern_c;
            end
            if type == 0
                y = zeros(obj.n);
            else
                y = sym('V', [obj.n obj.n]);
                for i = 1:obj.n
                    for j = 1:obj.n
                        y(i,j) = '0';
                    end
                end                
            end
            %%%
            % This part of the code is concerned only with the diagonal part
            %%%
            if type == 0
                Th = zeros(1,obj.n-1);
                for i = 1:(obj.n-1)
                    if rc == 0
                       j = obj.T_r(i);
                    elseif rc == 1
                       j = obj.T_c(i); 
                    end
                    if j > 0
                        Th(j) = h(i);
                    end
                end
                y(1,1) = Th(1);
                y(obj.n,obj.n) = -Th(obj.n-1);
                for i = 2:(obj.n-1)
                   y(i,i) = Th(i) - Th(i-1);
                end
            else
                Th = sym('Th', [obj.n-1,1]);
                for i =1:(obj.n-1)
                    Th(i) = sym('0');
                end
                for i = 1:(obj.n-1)
                    if rc == 0
                       j = obj.T_r(i);
                    elseif rc == 1
                       j = obj.T_c(i); 
                    end
                    if j > 0
                        Th(j) = h(i);
                    end
                end
                y(1,1) = Th(1);
                y(obj.n,obj.n) = -Th(obj.n-1);
                for i = 2:(obj.n-1)
                   y(i,i) = Th(i) - Th(i-1);
                end                
            end
            %%%
            % END of the part that deals with the diagonal
            %%%
            %%%
            % The following block of code just permutes the blocks given by
            % G1 and G2. It does not affect the diagonal, which has to be
            % treated separately 
            %%%
            i = 1;
            while true
                if rc == 0
                   j = obj.T_r(i);
                elseif rc == 1
                   j = obj.T_c(i); 
                end
               if j > 0
                   %calculate the size of the block that starts at (i,i)
                   %note: this code is adapted for the oriented case
                   %for a non-oriented BD triple, need a different one
                   size = 2;
                   for k=1:(obj.n-1-i)
                      if G1_pattern(i,i+k+1) ~= 0
                          size = size + 1;
                      else
                          break;
                      end
                   end
                   %now use the size to swap the blocks
                   for k=1:size
                       for m=1:size
                           if k~= m %to avoid changing the diagonal at this point
                               y(j+k-1,j+m-1) = x(i+k-1,i+m-1);
                           end
                       end
                   end
                   i = i + size;
                   if i > obj.n-2
                       break;
                   end
               elseif i > obj.n-2
                   break;
               else
                   i = i +1;
                   continue;
               end
            end
            %%%
            % End of pertmutation block
            %%%
        end

        function y = Tinv(obj,x,rc)
            type = strcmp(class(x),'sym');
            x = poisson_bracket.proj_sln(x);
            h = obj.cartan_coeff(x);
            if rc == 0
               G2_pattern = obj.G2_pattern_r;
            elseif rc == 1
               G2_pattern = obj.G2_pattern_c;
            end
            if type == 0
                y = zeros(obj.n);
            else
                y = sym('V', [obj.n obj.n]);
                for i = 1:obj.n
                    for j = 1:obj.n
                        y(i,j) = '0';
                    end
                end               
            end
            %%%
            % This part of the code is concerned only with the diagonal part
            %%%
            if type == 0
                Th = zeros(1,obj.n-1);
                for i = 1:(obj.n-1)
                    if rc == 0
                        j = obj.Tinv_r(i);
                    elseif rc == 1
                        j = obj.Tinv_c(i);
                    end
                    if j > 0
                        Th(j) = h(i);
                    end
                end
                y(1,1) = Th(1);
                y(obj.n,obj.n) = -Th(obj.n-1);
                for i = 2:(obj.n-1)
                   y(i,i) = Th(i) - Th(i-1);
                end
            else
                Th = sym('Th', [obj.n-1,1]);
                for i =1:(obj.n-1)
                    Th(i) = sym('0');
                end
                for i = 1:(obj.n-1)
                    if rc == 0
                        j = obj.Tinv_r(i);
                    elseif rc == 1
                        j = obj.Tinv_c(i);
                    end
                    if j > 0
                        Th(j) = h(i);
                    end
                end
                y(1,1) = Th(1);
                y(obj.n,obj.n) = -Th(obj.n-1);
                for i = 2:(obj.n-1)
                   y(i,i) = Th(i) - Th(i-1);
                end                
            end
            %%%
            % END of the part that deals with the diagonal
            %%%
            %%%
            % The following block of code just permutes the blocks given by
            % G1 and G2. It does not affect the diagonal, which has to be
            % treated separately 
            %%%
            i = 1;
            while true
                if rc == 0
                    j = obj.Tinv_r(i);
                elseif rc == 1
                    j = obj.Tinv_c(i);
                end
               if j > 0
                   %calculate the size of the block that starts at (i,i)
                   %note: this code is adapted for the oriented case
                   %for a non-oriented BD triple, need a different one
                   size = 2;
                   for k=1:(obj.n-1-i)
                      if G2_pattern(i,i+k+1) ~= 0
                          size = size + 1;
                      else
                          break;
                      end
                   end
                   %now use the size to swap the blocks
                   for k=1:size
                       for m=1:size
                           if k~= m %to avoid changing the diagonal at this point
                               y(j+k-1,j+m-1) = x(i+k-1,i+m-1);
                           end
                       end
                   end
                   i = i + size;
                   if i > obj.n-2
                       break;
                   end
               elseif i > obj.n-2
                   break;
               else
                   i = i + 1;
                   continue;
               end
            end
            %%%
            % End of pertmutation block
            %%%
        end

        function y = TTr(obj,x)
           y = obj.TT(x,0); 
        end
        function y = TTc(obj,x)
           y = obj.TT(x,1); 
        end
        function y = TTinvr(obj,x)
           y = obj.TTinv(x,0); 
        end
        function y = TTinvc(obj,x)
           y = obj.TTinv(x,1); 
        end
        
        function y = TT(obj,x,rc)
            %realizes the function (1-T)^{-1}
            %y = x;%zeros(obj.n);
            y = poisson_bracket.proj_sln(x);
            type = strcmp(class(x),'sym');
            %это вообще зачем было? зачем единицы на диагонале? 
%             for i =1:obj.n
%                y(i,i) = 1; 
%             end
            if type == 0
                zero = zeros(obj.n);
            else
                zero = sym('z',[obj.n obj.n]);
                for i = 1:obj.n
                    for j = 1:obj.n
                        zero(i,j) = sym('0');
                    end
                end
            end
            t = 1;
            A = obj.T(x,rc);
            while true
               y = y + A;
               A = obj.T(A,rc);
               t = t + 1;
               if (A == zero)
                   break;
               end
               if (t > 10000)
                  fprintf('The number of iterations exceeded.\n');
                  fprintf('Review the code for computing (1-T)^{-1}.\n');
                  break;
               end
            end
        end

        function y = TTinv(obj,x,rc)
            %realizes the function (1-T)^{-1}
            y = poisson_bracket.proj_sln(x);%zeros(obj.n);
            type = strcmp(class(x),'sym');
            if type == 0
                zero = zeros(obj.n);
            else
                zero = sym('z',[obj.n obj.n]);
                for i = 1:obj.n
                    for j = 1:obj.n
                        zero(i,j) = sym('0');
                    end
                end                
            end
            t = 1;
            A = obj.Tinv(x,rc);
            while true
               y = y + A;
               A = obj.Tinv(A,rc);
               t = t + 1;
               if (A == zero)
                   break;
               end
               if (t > 10000)
                  fprintf('The number of iterations exceeded.\n');
                  fprintf('Review the code for computing (1-T)^{-1}.\n');
                  break;
               end
            end
        end
    end
    methods (Static)
       function t = iselt(k,A)
            t = 0;
            for i = 1:length(A)
               if A(i) == k
                   t = 1;
                   return;
               end
            end
        end
       function G_pattern = ini_G_pattern(G,n)
            G_pattern = zeros(n,n);
            for i = 1:(n-1)
                h = extended_T.iselt(i,G);
                if h == 0
                    continue;
                else
                    G_pattern(i,i) = 1;
                    G_pattern(i+1,i+1) = 1;
                    G_pattern(i,i+1) = 1;
                    G_pattern(i+1,i) = 1;
                end
                for j = (i+1):(n-1)
                    h = extended_T.iselt(j,G);
                    if h == 1
                        for k=i:(j+1)
                            G_pattern(k,j+1) = 1;
                            G_pattern(j+1,k) = 1;
                        end
                    elseif h == 0
                       break;
                    end
                end
            end
        end
    end
end

