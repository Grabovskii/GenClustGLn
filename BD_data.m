%C:\Users\ДмитрийАлександрович\Desktop\УЧЕБА\еще учебники\Cluster Algebras\Plethora_generator

classdef BD_data < handle
    %this class represents a BD-triple and basic methods to work with it
    
    properties
        G1_r, G2_r, n %
        G1_c, G2_c,
        len_r %length of G_r
        len_c %length of G_c
    end
    
    methods
        function obj = BD_data(G1_r,G2_r,G1_c,G2_c,n)
            if nargin == 5
                obj.G1_r = G1_r; 
                obj.G2_r = G2_r;
                obj.G1_c = G1_c; 
                obj.G2_c = G2_c;
                obj.n = n;
                if G1_r(1) ~= 0
                    obj.len_r = length(G1_r);
                else
                   obj.len_r = 0; 
                end
                if G1_c(1) ~= 0
                    obj.len_c = length(G1_c);
                else
                   obj.len_c = 0; 
                end
            elseif nargin == 3
                obj.G1_r = G1_r;
                obj.G2_r = G2_r;
                obj.G1_c = obj.G1_r;
                obj.G2_c = obj.G2_r;
                obj.n = G1_c; %yes, G1_c is n if nargin == 3
                if G1_r(1) ~= 0
                    obj.len_r = length(G1_r);
                    obj.len_c = obj.len_r;
                else
                   obj.len_r = 0; 
                   obj.len_c = 0;
                end
            end
        end
        
        function r = T_(obj,i,rc)
            if strcmpi(rc,'r')
                r = obj.T_r(i);
            elseif strcmpi(rc,'c')
                r = obj.T_c(i);
            end
        end
        
        function r = Tinv_(obj,i,rc)
            if strcmpi(rc,'r')
                r = obj.Tinv_r(i);
            elseif strcmpi(rc,'c')
                r = obj.Tinv_c(i);
            end
        end
        
        function r = T_r(obj,i)
            r = -1;
            for k = 1:obj.len_r
               if (obj.G1_r(k) == i)
                  r = obj.G2_r(k);
                  return;
               end
            end
        end
        function r = Tinv_r(obj,i)
            r = -1;
            for k = 1:obj.len_r
               if obj.G2_r(k) == i
                   r = obj.G1_r(k);
                   return;
               end
            end
        end
        function r = T_c(obj,i)
            r = -1;
            for k = 1:obj.len_c
               if (obj.G1_c(k) == i)
                  r = obj.G2_c(k);
                  return;
               end
            end
        end
        function r = Tinv_c(obj,i)
            r = -1;
            for k = 1:obj.len_c
               if obj.G2_c(k) == i
                   r = obj.G1_c(k);
                   return;
               end
            end
        end
        function i_n = in(obj,G,i)
            n = obj.n;
            l = length(G);
            G_sorted = sort(G);

            i_n = 0;
            k = 1;

            for q=1:(i-1)
                if k <= l
                    if (G_sorted(k) == q)
                       k = k + 1; 
                    else
                       i_n = q;
                    end 
                else
                    i_n = q;
                end
            end
        end
        function i_p = ip(obj,G,i)
            %n = obj.n;
            l = length(G);
            G_sorted = sort(G);
            i_p = obj.n;
            k = l;
            for q = 1:(obj.n-i)
                if k > 0
                    if (G_sorted(k) == obj.n-q)
                        k = k - 1;
                    else
                       i_p = obj.n-q; 
                    end
                else 
                   i_p = obj.n-q; 
                end
            end
        end
        function [i_n,i_p] = ipm(obj,G,i)
            G_sorted = sort(G);
            i_n = obj.in(G_sorted,i);
            i_p = obj.ip(G_sorted,i);
            %This function returns the i minus and i plus values of i
            %to get the corresponding run, one need to add one to i_n
        end
        function output = runs(obj,G)
            %UNTITLED5 Summary of this function goes here
            %   Detailed explanation goes here
            output = zeros(obj.n,2);
            i=1;
            k = 1;
            while k <= obj.n
                [x,y] = obj.ipm(G,k);
                output(i,1) = x+1;
                output(i,2) = y;
                k = y + 1;
                i = i + 1;
            end
            output = output(1:(i-1),1:2);
        end
        
        function output = druns(obj,rc)
            %double runs: a run in G1 and the corresponding run in G2
            if strcmpi(rc,'c')
               G1 = obj.G1_c;
               G2 = obj.G2_c;
            elseif strcmpi(rc,'r')
               G1 = obj.G1_r;
               G2 = obj.G2_r;
            end
            output = zeros(obj.n,4);
            i=1;
            k = 1;
            while k < obj.n
                [x,y] = obj.ipm(G1,k);
                if y-x == 1 %we won't remember for this task the trivial runs
                    k = k +1;
                    continue;
                end
                output(i,1) = x+1;
                output(i,2) = y;
                [x_,y_] = obj.ipm(G2,obj.T_(k,rc));
                output(i,3) = x_+1;
                output(i,4) = y_;
                i = i + 1;
                k = y + 1;
            end
            output = output(1:(i-1),1:4);
        end
        
        function cc = max_paths(obj)
            n = obj.n;
            cc = zeros(2*n,2*n+2); %the precise upper bound for length of paths is 2*n
            num = 1; %the number of the current max alt path
            permit = ones(2,n-1); %if permit(i,j) = 1, can start a max alt path in ith half at j
            %determine where is not an end of a max'l alt path:
            for i = 1:n-1
               permit(1,i) = ~BD_data.iselt(i,obj.G1_c);
               permit(2,i) = ~BD_data.iselt(i,obj.G2_r);
            end
            for j = 1:2
                r = BD_data.nonzero(permit(j,:));
                while r > 0
                   path = obj.path(r,j);
                   if isempty(path) %if we have cycles
                       cc = cell(0);
                       return;
                   end
                   permit(j,r) = 0;
                   len = length(path);
                   if j == 1
                       cc(num, 1:len) = path;
                   else
                       cc(num,1) = -1; cc(num,2) = -1;
                       cc(num,3:len+2) = path;
                   end
                   num = num + 1;
                   r = BD_data.nonzero(permit(j,:));
                end
            end
        end
        
        function p = path(obj,i,s)
           %s can be 1 or 2, which says in which half we start the path
           %this function returns an alternating maximal path that starts
           %at i
           path_w_loop = 1;
           n = obj.n;
           p = zeros(1,2*n-2);
           tally = ones(2,n-1); %to check that we have no cycles
           p(1) = i; 
           p(2) = n-i;
           tally(s,i) = 0;
           tally(s,n-i) = 0;
           k = 3; %current vertex
           if s == 1
               next = obj.T_r(n-i);
           elseif s == 2
               next = obj.Tinv_c(n-i);
           end
           while (next > 0)
               s = 2 - s + 1;
               if ~tally(s,next)
                  p = cell(0);
                  return;
               end
               p(k) = next;
               p(k+1) = n - next;
               tally(s,next) = 0;
               tally(s,n-next) = 0;
               if (next == n - next) && path_w_loop
                   path_w_loop = 0;
                   tally = ones(2,n-1);
               elseif (next == n - next) && ~path_w_loop
                   p = cell(0);
                   return;
               end
               if s == 1
                   next = obj.T_r(n-next);
               elseif s == 2
                   next = obj.Tinv_c(n-next);
               end
               k = k + 2;
           end
        end
    end
    methods(Static)
        function t = iselt(k,A)
            t = 0;
            for i = 1:length(A)
               if A(i) == k
                   t = 1;
                   return;
               end
            end
        end 
        function r = nonzero(A)
            %returns index of a first non-zero number in array A; if haven't found, -1
            r = -1;
            for k=1:length(A)
               if A(k) ~= 0
                    r = k;
                    return;
               end
            end
        end
    end
end

