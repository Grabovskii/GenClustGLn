classdef double_gln < plethora & BD_data & handle 

    properties
        G_PATTERN, H_PATTERN, P_PATTERN, F_PATTERN
        plethora_adj_matrix
    end
    properties (Access = private)
        NodeLabels, NodeColors, markers
        force_frozen %a cell array of variables that are frozen by the user
        num_force_frozen %how many has a user frozen
    end
    properties (Constant)
        %MATLAB_RELEASE = 'R2020b'; %this primarily affects node labels
        %MATLAB_RELEASE = 'R2017b';
    end
    
    methods
        function obj = double_gln(G1_r,G2_r,G1_c,G2_c,n)
            if nargin == 3
               n = G1_c;
               G1_c = G1_r;
               G2_c = G2_r;
            end
            obj@plethora(G1_r,G2_r,G1_c,G2_c,n);
            obj.dim = 2*n^2 - n +1;
            patterns = obj.initialize_patterns();
            obj.G_PATTERN = patterns{1};
            obj.H_PATTERN = patterns{2};
            obj.P_PATTERN = patterns{3};
            obj.F_PATTERN = patterns{4};
            %since I call plethora's constructor, I already have 
            %plethora's adjacency matrix:
            obj.plethora_adj_matrix = obj.adj_matrix;%obj.ini_adj();
            %redefine the adjacency matrix:
            obj.adj_matrix = obj.ini_double_adj(obj.plethora_adj_matrix);
            XYData = obj.initialize_coords();
            obj.XData = XYData(:,1);
            obj.YData = XYData(:,2);
            %initialize NodeColors,markers,NodeLabels for future use
            obj.NodeColors = zeros(obj.dim,3);
            obj.markers = cell(obj.dim,1);
            obj.num_frozen = 0;
            for i =1:obj.dim
                if obj.isfrozen(i)
                    obj.num_frozen = obj.num_frozen + 1;
                    obj.markers{i} = 's';
                    obj.NodeColors(i,:) = [0 0.4470 0.7410];%light blue: [0.3010 0.7450 0.9330];
                else
                    obj.markers{i} = 'o';
                    obj.NodeColors(i,:) = [0 0 0];
                end
            end
            obj.markers{obj.inp(1,1)} = 'p'; %special vertex
            obj.NodeColors(obj.inp(1,1),:) = [0.9290 0.6940 0.1250];
            obj.NodeLabels = obj.ini_node_labels();
            obj.force_frozen = {};
            obj.num_force_frozen = 0;
        end
        function mutate(obj,type,i,j)
            ind = obj.ind(i,j,type);
            D = zeros(obj.dim);
                for i = 1:obj.dim
                    for j = 1:obj.dim
                        if (obj.isfrozen(i))&&(obj.isfrozen(j)) %does it work?
                           continue; 
                        end
                        if i == ind || j == ind
                            D(i,j) = -obj.adj_matrix(i,j);
                        else
                            bij = obj.adj_matrix(i,j);
                            bik =  obj.adj_matrix(i,ind);
                            bkj =  obj.adj_matrix(ind,j);
                            D(i,j) = bij + (1/2)*(abs(bik)*bkj + bik*abs(bkj));
                        end
                    end
                end
                obj.adj_matrix = D;
        end

        function mutate_B(obj,s)
            n = obj.n;
            for i=n:-1:s
               obj.mutate('h',s,i)
            end
            for i = 1:n-s
               obj.mutate('f',i,n-s-i+1)
            end
            for i = s:-1:2
                obj.mutate('g',s,i);
            end
        end
        
        function freeze(obj,type,i,j)
           %freezes the variable. need only to remove the edges between frozen vars
           k = obj.num_force_frozen + 1;
           obj.num_force_frozen = k;
           new_force_frozen = cell(k,1);
           for u = 1:k-1
              new_force_frozen(u,1:3) = obj.force_frozen(u,1:3); 
           end
           obj.force_frozen(k,1:3) = {type,i,j};
           ii = obj.ind(i,j,type);
           %remove the arrows between frozen vertices
           for u = 1:obj.dim
              if obj.isfrozen(u)
                 obj.adj_matrix(ii,u) = 0; 
                 obj.adj_matrix(u,ii) = 0; 
              end
           end
           %update the color and the shape on the quiver itself
           obj.markers{ii} = 's';
           obj.NodeColors(ii,:) = [0 0.4470 0.7410];%light blue: [0.3010 0.7450 0.9330];
        end
        function b = isfrozen(obj,i)
            t = obj.back_d(i);
            b = 0;
            if strcmp(t{1},'g')||strcmp(t{1},'h')
                b = obj.is_frozen(obj.ind_q(t{2},t{3}),obj.plethora_adj_matrix);
            end
            if ~b
                for u=1:obj.num_force_frozen
                    if strcmp(obj.force_frozen{u,1},t{1})&&(t{2} == obj.force_frozen{u,2})&&(t{3} == obj.force_frozen{u,3})
                       b = 1;
                       return;
                    end
                end
            end
        end
        function plot(obj)
            %mat_release = strcmp('R2020b',obj.MATLAB_RELEASE);
            D_clear = quiv.clear(obj.adj_matrix);
            D_arced = zeros(obj.dim);
            D_straight = zeros(obj.dim);
            types = {'g','h','p','f'};
            %the following long loop just split D_clear
            %into D_arced and D_straight, so that we know where to draw
            %edges and where arcs
            for k1 = 1:4
                for k2 = 1:4
                    type1 = types{k1};
                    type2 = types{k2};
                    for i1 = 1:obj.n
                        for j1 = 1:obj.n
                            for i2 = 1:obj.n
                                for j2 = 1:obj.n      
                            if obj.permit(type1,i1,j1)&&obj.permit(type2,i2,j2)
                                ind1 = obj.ind(i1,j1,type1);
                                ind2 = obj.ind(i2,j2,type2);
                                dist = obj.metric(type1,i1,j1,type2,i2,j2);
                                if dist > 2
                                	D_arced(ind1,ind2) = D_clear(ind1,ind2);
                                elseif dist > 0
                                   D_straight(ind1,ind2) = D_clear(ind1,ind2);
                                end
                            end
                                end
                            end
                        end
                    end
                end
            end
            %let's draw
            G = digraph(D_straight);
            %line for R2017 release:
            %obj.NodeLabels = ones(1,length(obj.NodeLabels));
            %if mat_release
            obj.p = plot(G, 'XData', obj.XData, 'YData', obj.YData, 'Marker', obj.markers, 'MarkerSize',10,'EdgeColor','black', 'NodeColor', obj.NodeColors, 'NodeLabel',obj.NodeLabels,'NodeFontSize',13);
            %else
            %obj.p = plot(G, 'XData', obj.XData, 'YData', obj.YData, 'Marker', obj.markers, 'MarkerSize',10,'EdgeColor','black', 'NodeColor', obj.NodeColors, 'NodeLabel',obj.NodeLabels);
            %end
             
            
             hold on;
             obj.arcs(D_arced);
             hold off;
             %some extra decorations
             for l=2:(obj.n-2)
                 if abs(obj.adj_matrix(obj.inp(1,l),obj.inp(l+1,1))) > 0
                    highlight(obj.p,obj.inp(1,l),obj.inp(l+1,1),'LineStyle','--');
                 end
             end
             for i =1:obj.dim
                 for j =1:obj.dim
                    if obj.adj_matrix(i,j) > 1
                         highlight(obj.p,i,j,'LineWidth',1.5);
                    end
                    if obj.adj_matrix(i,j) > 2
                         highlight(obj.p,i,j,'LineWidth',2);
                    end
                end
             end
        end
        
        function D = ini_double_adj(obj,M)
            %given the adjacency matrix for Plethora's quiver,
            %this function produces the adjacency matrix for the quiver in the
            %double
            %the PATTERNS should be already initialized before calling this
            %function
            %M = obj.ini_adj();
            
            n = obj.n;
            N = obj.dim;

            D = zeros(N,N);

            %инициализируем внутрибереговые связи по g
            for i =1:n
                for j=1:i
                    for u =1:n
                        for v=1:u
                            h = M(obj.ind_q(i,j),obj.ind_q(u,v));
                            %this chunk of if code is to take into account
                            %that arcs going by default into g(3,3) 
                            %have to go to the other shore into h(3,3)
                            if ((u == n)&&(v == n)) || ((i == n)&&(j==n))
                                if abs(u-i)+abs(v-j) < 3
                                    D(obj.ing(i,j), obj.ing(u,v)) = h;
                                else
                                    if (u == n)&&(v==n)
                                        D(obj.ing(i,j),obj.inh(u,v)) = h;
                                    elseif (i == n)&&(j == n)
                                        D(obj.inh(i,j),obj.ing(u,v)) = h;
                                    end
                                end
                            else
                                D(obj.ing(i,j), obj.ing(u,v)) = h;
                            end
                            %D(obj.ing(u,v), obj.ing(i,j)) = -h;
                        end
                    end
                end
            end

            %инициализируем внутрибереговые связи по h
            for i = 1:n
                for j = i:n
                    for u = 1:n
                        for v = u:n
                            if ~((i == j)&&(u==v))
                                h = M(obj.ind_q(i,j),obj.ind_q(u,v));
                                if ((u == n)&&(v == n)) || ((i == n)&&(j==n))
                                    if abs(u-i)+abs(v-j) < 3
                                        D(obj.inh(i,j), obj.inh(u,v)) = h;
                                    else
                                        if (u == n)&&(v==n)
                                            D(obj.inh(i,j),obj.ing(u,v)) = h;
                                        elseif (i == n)&&(j == n)
                                            D(obj.ing(i,j),obj.inh(u,v)) = h;
                                        end
                                    end
                                else
                                    D(obj.inh(i,j), obj.inh(u,v)) = h;
                                end
                                %D(obj.inh(u,v), obj.inh(i,j)) = -h;
                            end
                        end
                    end
                end
            end

            %инициализируем строго межбереговые связи (исключая, конечно, диагонали)
            for i = 1:n
                for j = 1:i %(i,j) из левого берега
                    for u = 1:n
                        for v = u:n %(u,v) из правого берега
                            h = M(obj.ind_q(i,j),obj.ind_q(u,v));
                            if (u~=v)&&(i~=j) %связи на диагоналях остаются на своих берегах
                                D(obj.ing(i,j), obj.inh(u,v)) = h;
                                D(obj.inh(u,v), obj.ing(i,j)) = -h;
                            end
                            %D(obj.inh(u,v), obj.ing(i,j)) = -h;
                        end
                    end
                end
            end

            %начинаем вставлять ромб

            %special vertex
            D(obj.inp(1,1),obj.inh(1,1)) = 1;
            D(obj.inp(1,1),obj.inp(2,1)) = 1;
            D(obj.inp(1,1),obj.ing(1,1)) = -1;
            D(obj.inp(1,1),obj.inp(1,2)) = -1;

        %     if n>=3
        %        D(obj.inp(1,2),obj.inp(1,1)) = 1;
        %        D(obj.inp(1,1),obj.inp(1,2)) = -1;
        %        D(obj.inp(1,1),obj.inp(2,1)) = 1;
        %        D(obj.inp(2,1),obj.inp(1,1)) = -1;
        %     end

            %boundary case for phi
            for k=2:(n-2)
                l = n-k;
                D(obj.inp(k,l),obj.inp(k,l-1)) = 1;
                D(obj.inp(k,l), obj.inf(k,l-1)) = -1;
                D(obj.inp(k,l),obj.inf(k-1,l)) = 1;
                D(obj.inp(k,l),obj.inp(k-1,l)) = -1;
            end

            %case of f's
            for k = 1:(n-2)
               for l =1:(n-1-k)
                  D(obj.inf(k,l),obj.inf(k-1,l)) = 1;
                  D(obj.inf(k,l),obj.inf(k-1,l+1)) = -1;
                  D(obj.inf(k,l),obj.inf(k,l-1)) = -1;
                  D(obj.inf(k,l),obj.inf(k+1,l)) = -1;
                  D(obj.inf(k,l),obj.inf(k+1,l-1)) = 1;
                  D(obj.inf(k,l),obj.inf(k,l+1)) = 1;
               end
            end

            %case of g11, hnn, h11
            D(obj.ing(1,1),obj.inp(1,1)) = 1;
            D(obj.ing(1,1),obj.ing(2,2)) = 1;
            D(obj.ing(1,1),obj.inp(n-1,1)) = -1;

            D(obj.inh(n,n), obj.inf(1,1)) = -1;
            D(obj.inh(n,n),obj.ing(n,n)) = 1;
            D(obj.inh(n,n),obj.inh(n-1,n)) = 1;

            D(obj.inh(1,1),obj.inp(1,n-1)) = 1;
            D(obj.inh(1,1),obj.inp(1,1)) = -1;

            %case of hii for i~=1 and i~= n
            for i=2:(n-1)
                D(obj.inh(i,i), obj.inf(1,n-i)) = 1;
                D(obj.inh(i,i),obj.inh(i-1,i)) = 1;
                D(obj.inh(i,i),obj.inh(i,i+1)) = -1;
                D(obj.inh(i,i),obj.inf(1,n-i+1)) = -1;
            end

            %case of phi_{n-1, 1}
            D(obj.inp(n-1,1),obj.ing(2,2)) = -1;
            D(obj.inp(n-1,1),obj.ing(1,1)) = 1;
            D(obj.inp(n-1,1),obj.inp(1,n-2)) = -1;
            D(obj.inp(n-1,1),obj.inp(1,n-1)) = 1;
            D(obj.inp(n-1,1),obj.inf(n-2,1)) = 1;

            %case of phi_{1,n-1}
            D(obj.inp(1,n-1),obj.inp(1,n-2)) = 1;
            D(obj.inp(1,n-1),obj.inh(2,2)) = 1;
            D(obj.inp(1,n-1),obj.inf(1,n-2)) = -1;
            D(obj.inp(1,n-1),obj.inp(n-1,1)) = -1;
            D(obj.inp(1,n-1),obj.inh(1,1)) = -1;

            %phi_{k,1} for k in [2,n-2]
            for k = 2:(n-2)
               D(obj.inp(k,1),obj.inp(k-1,2)) = 1;
               D(obj.inp(k,1),obj.inp(k,2)) = -1;
               D(obj.inp(k,1),obj.inp(1,k-1)) = -1;
               D(obj.inp(k,1),obj.inp(1,k)) = 1;
            end
            %phi_{1,l} for l in [2,n-2]
            for l = 2:(n-2)
               D(obj.inp(1,l),obj.inp(l,1)) = -1;
               D(obj.inp(1,l),obj.inp(2,l)) = 1;
               D(obj.inp(1,l),obj.inp(1,l-1)) = 1;
               D(obj.inp(1,l),obj.inp(1,l+1)) = -1;
               D(obj.inp(1,l),obj.inp(l+1,1)) = 1;
               D(obj.inp(1,l),obj.inp(2,l-1)) = -1;
               D(obj.inp(1,l),obj.inp(l,1)) = -1;
            end

            %interior phi's, i.e. for k,l ~=1 and k+l < n
            for k=2:(n-2)
                for l=2:(n-1-k)
                    D(obj.inp(k,l),obj.inp(k,l+1)) = -1;
                    D(obj.inp(k,l),obj.inp(k+1,l)) = 1;
                    D(obj.inp(k,l),obj.inp(k+1,l-1)) = -1;
                    D(obj.inp(k,l),obj.inp(k,l-1)) = 1;
                    D(obj.inp(k,l),obj.inp(k-1,l+1)) = 1;
                    D(obj.inp(k,l),obj.inp(k-1,l)) = -1;
                end
            end

            %case g_{i,i} i ~= n, i~= 1
            for i = 2:(n-1)
                D(obj.ing(i,i),obj.ing(i,i+1)) = -1;
                D(obj.ing(i,i),obj.ing(i+1,i+1)) = 1;
                D(obj.ing(i,i),obj.ing(i+1,i)) = -1;
                D(obj.ing(i,i),obj.ing(i,i-1)) = 1;
                D(obj.ing(i,i),obj.ing(i-1,i-1)) = -1;
                D(obj.ing(i,i),obj.ing(i-1,i)) = 1;
            end

            %case g_{n,n}
            D(obj.ing(n,n),obj.ing(n,n+1)) = -1;
            D(obj.ing(n,n),obj.ing(n,n-1)) = 1;
            D(obj.ing(n,n),obj.ing(n-1,n-1)) = -1;
            D(obj.ing(n,n),obj.ing(n-1,n)) = 1;

            %handle separately the case with the double edge
            if n > 3
               D(obj.inp(2,1),obj.inp(1,2)) = 2; 
               D(obj.inp(1,2),obj.inp(2,1)) = -2;
            end

        end
        function k = ind(obj,i,j,s)
            %calculates the index in the adjacency matrix
            %using indices of g, h, f, and phi parts
            if strcmp(s,'g')
                k = obj.G_PATTERN(i,j);
            elseif strcmp(s,'h')
                k = obj.H_PATTERN(i,j);
            elseif strcmp(s,'f')
                k = obj.F_PATTERN(i,j);
            elseif strcmp(s,'p')||strcmp(s,'phi')
                k = obj.P_PATTERN(i,j);
            end
        end
        function k = ing(obj,i,j)
            if j == i+1
                k = obj.inf(obj.n-i,1);
            else
                k = obj.ind(i,j,'g');
            end
        end
        function k = inh(obj,i,j)
            k = obj.ind(i,j,'h');
        end
        function k = inp(obj,i,j)
            k = obj.ind(i,j,'p');
        end
        function k = inf(obj,i,j)
            n = obj.n;
            %we realize here Rermark 3.1 from page 8
            if (i+j == n)&&(i > 0)&&(j > 0)
               k = obj.ind(i,j,'p'); 
            elseif (i == 0)
               k = obj.inh(n-j+1,n-j+1); 
            elseif (j == 0)
               k = obj.ing(n-i+1,n-i+1); 
            else
               k = obj.ind(i,j,'f');
            end
        end 
        function Data = initialize_coords(obj)
            n = obj.n;
            N = obj.dim;
            XData = zeros(N,1);
            YData = zeros(N,1);

            %assign coordinates to g's shore
            for i =1:n
                for j = 1:i
                   YData(obj.ing(i,j)) = (n-i+1); %заметаем от 1 до n, т.е. фи пойдет с n+1
                   XData(obj.ing(i,j)) = (n-i)*(1/2) + j;
                end
            end

            %assign coordinates to h's shore
            for i =1:n
                for j =i:n
                    YData(obj.inh(i,j)) = (n-j+1);
                    %if read the picture from left to right, get
                    XData(obj.inh(i,j)) = 1 + (n-1) + 1 + (j-i)*(1/2) + (n-i)*(1/2);
                end
            end

            %assign coordinates to phi's part of the rhomb
            for i =1:(n-1)
                for j=1:(n-i)
                    %here (i+j) is constant on horizontal slices ...
                    %reads: 1 + n + (n-1-i)
                    YData(obj.inp(i,j)) = 2*n-2 - (i+j-2);
                    XData(obj.inp(i,j)) = n+(1/2) - (1/2)*(i-j);%nn-(1/2) - (1/2)*((i-j) - (nn-2));
                end
            end
            %assign coordinates to f's part of the rhomb
            for i =1:(n-2)
                for j=1:(n-i-1)
                    YData(obj.inf(i,j)) = i+j;
                    XData(obj.inf(i,j)) = n+(1/2) - (1/2)*(i-j); %nn-(1/2)*(i-j) + (nn-3)/2;
                end
            end

            Data = zeros(N,2);
            Data(:,1) = XData;
            Data(:,2) = YData;
        end      
        function NodeLabels = ini_node_labels(obj)
            %mat_release = strcmp('R2020b',obj.MATLAB_RELEASE);
            %if mat_release 
            b = '_'; 
            phi = '\phi';
            %else
            %    b = '';
            %    phi = 'p';
            %end
            %latex or tex?
            %dol = '$';
            dol = '';
            n = obj.n;
            N = obj.dim;
            NodeLabels = cell(N,1);
            for i = 1:n
               for j = 1:i
                   NodeLabels{obj.ing(i,j)} = append(dol,'g',b, '{',num2str(i), num2str(j),'}',dol);
               end
            end
            for i = 1:n
                for j = i:n
                    NodeLabels{obj.inh(i,j)} = append(dol,'h',b, '{', num2str(i), num2str(j),'}',dol);
                end
            end
            for i = 1:(n-2)
                for j = 1:(n-i-1)
                    NodeLabels{obj.inf(i,j)} = append(dol,'f',b, '{',num2str(i), num2str(j),'}',dol);
                end
            end
            for i = 1:(n-1)
                for j = 1:(n-i)
                    NodeLabels{obj.inp(i,j)} = append(dol,phi,b,'{', num2str(i), num2str(j),'}',dol);
                end
            end
        end
        function results = initialize_patterns(obj)
            results = cell(4);
            n = obj.n;
            results{1} = zeros(n,n);%G_PATTERN
            results{2} = zeros(n,n); %H_PATTERN
            results{3} = zeros(n,n); %P_PATTERN
            results{4} = zeros(n,n); %F_PATTERN

            %заметает нижний левый угол матрицы
            s = 0;
            for i =1:n
                for j = 1:i
                    s = s + 1;
                    results{1}(i,j) = s;
                end
            end

            %заметает правый верхний угол матрицы
            for i = 1:n
                for j =i:n
                    s = s+1;
                    results{2}(i,j) = s;
                end
            end

            %левый верхний угол, строго
            for i=1:(n-1) %here omit the diagonal
                for j = 1:n-i
                    s = s+1;
                    results{3}(i,j) = s;
                end
            end

            %левый верхний, дважды строго
            for i =1:(n-2)
                for j = 1: n-i - 1
                    s = s+1;
                    results{4}(i,j) = s;
                end
            end
        end
        function b = back_d(obj,ind)
           %the function eats the index of a function in terms of the
           %indexation in the adjoint matrix and returns
           %a cell array b = {'type',i,j},
           %where 'type' = 'g' or 'h' or 'f' or 'phi'
           %and (i,j) is the index of the function
           b = cell(3);
           for i = 1:obj.n
               for j = 1:obj.n
                    if obj.G_PATTERN(i,j) == ind
                        b{1} = 'g';
                        b{2} = i;
                        b{3} = j;
                        return;
                    elseif obj.H_PATTERN(i,j) == ind
                        b{1} = 'h';
                        b{2} = i;
                        b{3} = j;
                    elseif obj.F_PATTERN(i,j) == ind
                        b{1} = 'f';
                        b{2} = i;
                        b{3} = j;
                    elseif obj.P_PATTERN(i,j) == ind
                        b{1} = 'phi';
                        b{2} = i;
                        b{3} = j;
                    end
               end
           end
        end
    end
    
    methods(Access = private)
        function p = permit(obj,type,i,j)
           p = strcmp(type,'g') && j <= i;
           p = p || (strcmp(type,'h') && j >= i);
           p = p || (strcmp(type,'f') && (i+j <= obj.n-1));
           p = p || ((strcmp(type,'p') || strcmp(type,'phi'))&& i+j <= obj.n);
        end
        function d = metric(obj,type1,i1,j1,type2,i2,j2)
           d = -1;
           %WHEN type1 and type2 ARE THE SAME
           if strcmp(type1,type2)
               % I want an exception for phi's on the sides
               % so that I draw dashed edges between them
               if strcmp(type1,'g')&& (abs(i1-i2)==2) && (abs(j1-j2)==1)
                   d = 1;
                   return;
               end
               if strcmp(type1,'phi')||strcmp(type1,'p')%||strcmp(type1,'f')
                   c = (i1 == 1 && j2 == 1)||(j1 == 1 && i2 == 1);
                   c = c&&((abs(i1-j2) == 1)||(abs(i2-j1) == 1));
                    if c
                        d = 1;
                        return;
                    end
               end
               if (abs(i1-i2) == 1) && (abs(j1-j2) == 2)
                  d = 1;
                  return;
               end
               d = abs(i1-i2) + abs(j1-j2);
               if j1==j2
                   d = d+1;
               end
               return;
           %WHEN type1 and type2 ARE DIFFERENT
           elseif strcmp(type1,'g') && strcmp(type2,'h')
              if i1 == i2 && j1 == j2 && i1 == obj.n && i2 == obj.n
                  d = 1;
              elseif i1 == obj.n && j1 == obj.n && i2 == obj.n-1 && j2 == obj.n-1
                  d =1;
              else
                  d = 3;
              end
           elseif strcmp(type1,'g') && strcmp(type2, 'f')
              if ((i1+i2 == obj.n)&&(j2==2)) || ((i1+i2 == obj.n+1)&&(j2 ==1))
                  d = 1;
              elseif (i1 ~= j1)||(j2 ~= 1)
                  d = 3;
              elseif abs(i1-(obj.n-i2)) < 2
                  d = 1;
              else
                  d = 3;
              end
           elseif strcmp(type1,'h') && strcmp(type2,'f')
              if (i1 ~= j1)||(i2 ~= 1)
                  d =3;
              elseif abs(j1-(obj.n-j2)) < 2
                  d = 1;
              else
                  d= 3;
              end
           elseif strcmp(type1,'g') && (strcmp(type2,'phi') || strcmp(type2,'p'))
              if (i1 == j1)&&(i1 == 2)&&(i2 + j2 == obj.n)
                  d = 1;
              elseif i1== 1 && j1 == 1 && j2 == 1
                  d = 1;
              else 
                  d = 3;
              end
           elseif (strcmp(type1,'p') || strcmp(type1,'phi'))&&(strcmp(type2,'f'))
               if (i1+j1 == obj.n)&&(i2+j2 == obj.n-1)
                   if (abs(j1-j2) + abs(i1-i2) <= 3)
                      d = 1; 
                   else
                      d = abs(i1-i2) + abs(j1-j2);
                   end
               else
                   d = 3;
               end
           elseif (strcmp(type1,'h'))&&((strcmp(type2,'p'))||(strcmp(type2,'phi')))
              if (i1==j1) && (i1 < 3) && (i2 == 1) && (j2 == obj.n-1)
                  d = 1;
              elseif i1== 1 && j1 == 1 && i2 ==1 && j2 == 1
                  d = 1;
              elseif (i1 == j1)&&(i1 == 2)&&(j2 == obj.n-2) && (i2 == 2)
                  d = 1;
              elseif (i1 == j1-1)&&(i1 == 1)&&(j2 == obj.n-1) && (i2 == 1)
                  d = 1;
              else
                  d = 3;
              end
           else% d < 0
              d = obj.metric(type2,i2,j2,type1,i1,j1); 
              %return;
           end
        end
    end 
end

