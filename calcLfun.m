classdef calcLfun < handle & BD_data
    %just a container
    
    %just some thoughts: I have to drag along all of this data (properties)
    %to avoid the following two alternatives:
    %either each calcLfun has a shorestech field, or I realize in calcL
    %all methods I need for functions, which also seems awkward
    %for then calcLfun.m drags weird field along
    %I would still need a receptacle for functions....
    
    properties
        type
        i
        j
        s
        fun
        Lm
        grL
        LgrL
        grLL
        
        L
        bL
        K
        bK
        Phi
        Psi
        
        b
        a
        bb
        ba
        p
        q
        nblocks
        exitX
        exitY
    end
    
    methods
        function obj = calcLfun(G1_r,G2_r,G1_c,G2_c,n)
           obj@BD_data(G1_r,G2_r,G1_c,G2_c,n);
        end
        
        
        function update(obj,varargin)
            for k = 1:length(varargin)
               if strcmpi(varargin{k},'type')
                  obj.type = varargin{k+1}; 
               elseif strcmpi(varargin{k},'i')
                   obj.i = varargin{k+1};
               elseif strcmpi(varargin{k},'j')
                   obj.j = varargin{k+1};
               elseif strcmpi(varargin{k},'fun')
                   obj.fun = varargin{k+1}; 
               elseif strcmpi(varargin{k},'L')
                   obj.L = varargin{k+1};
               elseif strcmpi(varargin{k},'K')
                   obj.K = varargin{k+1};
               elseif strcmpi(varargin{k},'Lm')
                   obj.Lm = varargin{k+1};
               elseif strcmpi(varargin{k},'bK')
                   obj.bK = varargin{k+1};
               elseif strcmpi(varargin{k},'Phi')
                   obj.Phi = varargin{k+1};
               elseif strcmpi(varargin{k},'Psi')
                   obj.Psi = varargin{k+1};
               elseif strcmpi(varargin{k},'bL')
                   obj.bL = varargin{k+1};
               elseif strcmpi(varargin{k},'a')
                   obj.a = varargin{k+1};
               elseif strcmpi(varargin{k},'b')
                   obj.b = varargin{k+1};
               elseif strcmpi(varargin{k},'bb')
                   obj.bb = varargin{k+1};
               elseif strcmpi(varargin{k},'ba')
                   obj.ba = varargin{k+1};
               elseif strcmpi(varargin{k},'p')
                   obj.p = varargin{k+1};
               elseif strcmpi(varargin{k},'q')
                   obj.q = varargin{k+1}; 
               elseif strcmpi(varargin{k},'s')
                   obj.s = varargin{k+1};
               elseif strcmpi(varargin{k},'nblocks')
                   obj.nblocks = varargin{k+1}; 
               elseif strcmpi(varargin{k},'exitX')
                   obj.exitX = varargin{k+1}; 
               elseif strcmpi(varargin{k},'exitY')
                   obj.exitY = varargin{k+1}; 
               end
            end
            
            [N,~] = size(obj.Lm);
            
            d = det(obj.fun)*inv(obj.fun);
            obj.grL = sym(zeros(N));
            obj.grL(obj.s:N,obj.s:N) = d;
            
            obj.LgrL = obj.Lm*obj.grL;
            obj.grLL = obj.grL*obj.Lm;
            
        end
        
        function result = geti(obj,grtype,interval1,varargin)
            %the same as get except allows arbitrary intervals
            %basically I need these functions to simplify the code
           if isempty(varargin)
               interval2 = interval1;
           else
               interval2 = varargin{1};
           end
           if ~isempty(interval1)&& ~isempty(interval2)
               %if interval1(1)*interval1(2) ~=0 && interval2(1)*interval2(2) ~= 0
               gr = obj.getGradient(grtype);
               result = gr(interval1(1):interval1(2),interval2(1):interval2(2));
               %end
           else
               result = 0;
           end
        end
        
        function result = get(obj,t,grtype,btype1,varargin)
            %rows first, then columns; so, btype1 for rows
            %if the input is in terms of types, then we use getInterval
           if isempty(varargin)
               btype2 = btype1;
           else
               btype2 = varargin{1}; 
           end
           ii1 = obj.getInterval(t,btype1);
           ii2 = obj.getInterval(t,btype2);
           if ~isempty(ii1) && ~isempty(ii2)
               gr = obj.getGradient(grtype);
               result = gr(ii1(1):ii1(2),ii2(1):ii2(2));
           else
              result = 0; 
           end
           
        end

        function r = getBound(obj,t,type)
           %here we get a,b,ba,bb,exitX,exitY
           r = 0;
           if strcmpi(type,'a')
               if (0 < t && t <= obj.nblocks)
                  r = obj.a(t); 
               end
           elseif strcmpi(type,'b')
               if (0 < t && t <= obj.nblocks+1)
                  r = obj.b(t); 
               end
           elseif strcmpi(type,'ba')
               if (0 < t && t <= obj.nblocks)
                  r = obj.ba(t); 
               end
           elseif strcmpi(type,'bb')
               if (0< t && t <= obj.nblocks)
                   r = obj.bb(t);
               end
           elseif strcmpi(type,'exitX')
               if (0< t && t <= obj.nblocks)
                   r = obj.exitX(t);
               end
           elseif strcmpi(type,'exitY')
               if (0< t && t <= obj.nblocks)
                   r = obj.exitY(t);
               end
           end
        end
        
        function ii = getInterval(obj,t,type)
           if strcmpi(type,'Psi')
               if (0 < t && t <= obj.nblocks+1)
                   ii = obj.Psi{t};
                   return;
               else
                   ii = cell(0); 
               end
           elseif strcmpi(type,'Psi+1')
               ii = obj.getInterval(t+1,'Psi');
               return;
           elseif strcmpi(type,'Phi')
               if (0 < t && t<=obj.nblocks)
                  ii = obj.Phi{t}; 
               else
                  ii = cell(0);
               end
           elseif strcmpi(type,'bK')
               if (0 < t && t <= obj.nblocks)
                   ii = obj.bK{t};
               else
                   ii = cell(0);
               end
           elseif strcmpi(type,'bK-1')
               ii = obj.getInterval(t-1,'bK');
           elseif strcmpi(type,'bL')
               if (0 <= t && t <= obj.nblocks)
                   %remember: there is a shift, for we want to have L{0}
                   ii = obj.bL{t+1};
               else
                   ii = cell(0);
               end
           elseif strcmpi(type,'K')
               if (0 < t && t <= obj.nblocks)
                   ii = obj.K{t};
               else
                   ii = cell(0);
               end
           elseif strcmpi(type,'L')
               if (0 < t && t <= obj.nblocks+1)
                   ii = obj.L{t};
               else
                   ii = cell(0);
               end
           elseif strcmpi(type,'L-Psi')
               L = obj.getInterval(t,'L');
               Psi = obj.getInterval(t,'Psi');
               if ~isempty(Psi) && ~isempty(L)
                  ii(1) = L(1);
                  ii(2) = Psi(1)-1;
                  if ii(1) > ii(2)
                      ii = cell(0);
                  end
               elseif isempty(Psi) && ~isempty(L)
                   ii = L;
               else
                  ii = cell(0); 
               end
           elseif strcmpi(type,'bL-Psi')
               bL = obj.getInterval(t,'bL');
               Psi = obj.getInterval(t+1,'Psi');
               if ~isempty(Psi) && ~isempty(bL)
                  ii(2) = bL(2);
                  ii(1) = Psi(2)+1;
               elseif ~isempty(bL) && isempty(Psi)
                  ii = bL;
               else
                  ii = cell(0); 
               end
           elseif strcmpi(type,'bK-Phi')
               bK = obj.getInterval(t,'bK');
               Phi = obj.getInterval(t,'Phi');
               if ~isempty(Phi) && ~isempty(bK)
                   ii(1) = bK(1);
                   ii(2) = Phi(1)-1;
                   if ii(1) > ii(2)
                       ii = cell(0); %this might happen if the intervals coincide
                   end
               elseif isempty(Phi) && ~isempty(bK)
                    ii = bK;
               else
                  ii = cell(0); 
               end
           elseif strcmpi(type,'K-Phi')
               K = obj.getInterval(t,'K');
               Phi = obj.getInterval(t,'Phi');
               if ~isempty(K) && ~isempty(Phi)
                   ii(2) = K(2);
                   ii(1) = Phi(2)+1;
                   if ii(1) > ii(2)
                       ii = cell(0);
                   end
               elseif isempty(Phi) && ~isempty(K)
                   ii = K;
               else
                   ii = cell(0);
               end
           end
       end
       function gr = getGradient(obj,type)
           if strcmpi(type,'grL')
               gr = obj.grL;
           elseif strcmpi(type,'LgrL')
               gr = obj.LgrL;
           elseif strcmpi(type,'grLL')
               gr = obj.grLL;
           elseif strcmpi(type,'L')
               gr = obj.Lm;
           end
       end
       
               
        function run = getRun(obj,t,type,runtype)
           %returns the run that corresponds to Phi or Psi of fun(num)
           %don't need for the main function, only for tests
           run = cell(0);
           if strcmpi(type,'Psi')
               bb = obj.getBound(t-1,'bb');
               b = obj.getBound(t,'b');
               if bb+b == 0
                   return;
               end
               if b ~= 0
                  run = zeros(1,2);
                  [run(1),run(2)] = obj.ipm(obj.G1_c, b);
                  run(1) = run(1)+1;
                  len = run(2)-run(1);
                  if strcmpi(runtype,'Y')
                      if len > 0
                        run(1) = obj.T_c(run(1));
                        run(2) = obj.T_c(run(2)-1)+1;
                      else
                        run = cell(0);
                      end
                  end
               elseif bb ~= 0
                  run = zeros(1,2);
                  [run(1),run(2)] = obj.ipm(obj.G2_c, bb);
                  run(1) = run(1)+1;
                  len = run(2)-run(1);
                  if strcmpi(runtype,'X')
                      if len > 0
                        run(1) = obj.Tinv_c(run(1));
                        run(2) = obj.Tinv_c(run(2)-1)+1;
                      else
                        run = cell(0);
                      end
                  end
               end
           end
           if strcmpi(type,'Phi')
               ba = obj.getBound(t,'ba');
               a = obj.getBound(t,'a');
               if ba+a == 0
                   return;
               end
               if a ~= 0
                  run = zeros(1,2);
                  [run(1),run(2)] = obj.ipm(obj.G1_r, a);
                  run(1) = run(1)+1;
                  len = run(2) - run(1);
                  if strcmpi(runtype,'Y')
                      if len > 0
                        run(1) = obj.T_r(run(1));
                        run(2) = obj.T_r(run(2)-1)+1;
                      else
                        run = cell(0);
                      end
                  end
               elseif ba ~= 0
                  run = zeros(1,2);
                  [run(1),run(2)] = obj.ipm(obj.G2_r, ba);
                  run(1) = run(1)+1;
                  len = run(2) - run(1);
                  if strcmpi(runtype,'X')
                      if len > 0
                          run(1) = obj.Tinv_r(run(1));
                          run(2) = obj.Tinv_r(run(2)-1)+1;
                      else 
                          run = cell(0);
                      end
                  end
               end
           end
        end
        
    end
end

