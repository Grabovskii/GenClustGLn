classdef calcL < handle & matrixClass & extended_T & genSettings
    %I create this class to relize methods for working with blocks
    %it requires switching from shores.m to shorestech.m in double_full.m
    %settings
    %need extended_T for projections onto G1,G2,etc.
    
    properties
        s %a copy of shores
        fun1
        fun2
    end
    
    methods
        function obj = calcL(shore)
            obj@extended_T(shore.G1_r,shore.G2_r,shore.G1_c,shore.G2_c,shore.n);
            obj.s = shore;
            obj.fun1 = calcLfun(shore.G1_r,shore.G2_r,shore.G1_c,shore.G2_c,shore.n);
            obj.fun2 = calcLfun(shore.G1_r,shore.G2_r,shore.G1_c,shore.G2_c,shore.n);
        end
        
        function set(obj,type,i,j,num)
            [Lind,s] = obj.s.Lindex(type,i,j);
            if strcmpi(type,'g')
                fun = obj.s.G{i,j};
            elseif strcmpi(type,'h')
                fun = obj.s.H{i,j};
            end
            a = obj.s.getii('a',Lind);
            b = obj.s.getii('b',Lind);
            ba = obj.s.getii('ba',Lind);
            bb = obj.s.getii('bb',Lind);
            nblocks = obj.s.getii('nblocks',Lind);
            p = obj.s.pIndex(type,i,j,'Lindex',Lind,'s',s);
            q = obj.s.qIndex(type,i,j,'Lindex',Lind,'s',s,'p',p);
            if num == 1
                obj.fun1.update('type',type,'i',i,'j',j,...
                    'Lm',obj.s.L{Lind},'L',obj.s.getii('L',Lind), ...
                    'K', obj.s.getii('K',Lind),'Phi',obj.s.getii('Phi',Lind),...
                    'Psi',obj.s.getii('Psi',Lind),'bL',obj.s.getii('bL',Lind),...
                    'bK',obj.s.getii('bK',Lind),'fun',fun,'b',b,'a',a,'p',p,...
                    'q',q,'bb',bb,'ba',ba,'s',s,'nblocks',nblocks,...
                    'exitX',obj.s.getii('exitX',Lind),'exitY',obj.s.getii('exitY',Lind));
            elseif num == 2
                obj.fun2.update('type',type,'i',i,'j',j,...
                    'Lm',obj.s.L{Lind},'L',obj.s.getii('L',Lind), ...
                    'K', obj.s.getii('K',Lind),'Phi',obj.s.getii('Phi',Lind),...
                    'Psi',obj.s.getii('Psi',Lind),'bL',obj.s.getii('bL',Lind),...
                    'bK',obj.s.getii('bK',Lind),'fun',fun,'b',b,'a',a,'p',p,...
                    'q',q,'bb',bb,'ba',ba,'s',s,'nblocks',nblocks,...
                    'exitX',obj.s.getii('exitX',Lind),'exitY',obj.s.getii('exitY',Lind));
            end
        end
        
        function r = sigma(obj,t,type,varargin)
            %blocks of the second functions are embedded into
            %the blocks of the first function
            %one letter for the second, two letters for the first
            
            %I need varargin to realize p-p-1 convention
            %for p-p convention, we only in Y_p;
            %for Plethora's p-p-1 convention, 
            %there could be either Y_p or Y_{p-1}
            r = cell(0);
            
            if ~isempty(varargin)
                index = varargin{1};
            else
                index = obj.fun1.p; %the index in p-p convention
                %in Plethora's article, they have obj.fun1.q, 
                %which is wrong
            end
            
           if strcmpi(type,'bK')
              KK = obj.fun1.getInterval(index,'bK');
              K = obj.fun2.getInterval(t,'bK');
              if ~isempty(K) && ~isempty(KK)
                r = [KK(1),KK(2)-(KK(2)-KK(1)+1)+(K(2)-K(1)+1)];
              else
                  r = cell(0);
                  return;
              end
           elseif strcmpi(type,'bL')
              KK = obj.fun1.getInterval(index,'bL');
              K = obj.fun2.getInterval(t,'bL');
              if ~isempty(K) && ~isempty(KK)
                r = [KK(1)+(KK(2)-KK(1)+1)-(K(2)-K(1)+1),KK(2)];
              else
                  r = cell(0);
                  return;
              end
           elseif strcmpi(type,'Phi')
              Phi = obj.fun2.getInterval(t,'Phi');
              if isempty(Phi)
                  r = cell(0);
                  return;
              end
              
%               if (index == obj.fun1.nblocks) && isempty(obj.fun1.getInterval(index,'bL'))
%                   %indicator of the empty last Y block;
%                   %then use the X-block to calculate the correct
%                   %embedding
%                   r = obj.rho(t,'K',index);
%                   if ~isempty(r)
%                        r(2) = r(1)+(Phi(2)-Phi(1));
%                        return;
%                   end
%               elseif (index == 1) && isempty(obj.fun1.getInterval(index,'L'))
%                   r = obj.sigma(t,'bK',index);
%                   if ~isempty(r)
%                       r(1) = r(2) - (Phi(2)-Phi(1));
%                       return;
%                   end
%               end
              
              r = obj.sigma(t,'bK',index);
              if ~isempty(r)
                  r(1) = r(2) - (Phi(2)-Phi(1));
              else
                  r = cell(0);
                  return;
              end
           elseif strcmpi(type,'Psi')
               r = obj.sigma(t-1,'bL',index);
               Psi = obj.fun2.getInterval(t,'Psi');
              if ~isempty(r) && ~isempty(Psi) 
                  r(2) = r(1)+(Psi(2)-Psi(1)); 
              else
                  r = cell(0);
              end
           elseif strcmpi(type,'bL-Psi')
              r = obj.sigma(t,'bL',index);
              Psi = obj.sigma(t+1,'Psi',index);
              if ~isempty(r) && isempty(Psi)
                  return;
              elseif ~isempty(r) && ~isempty(Psi)
                  r(1) = Psi(2)+1;
                  if r(1) > r(2)
                      r = cell(0);
                      return;
                  end
              else
                 r = cell(0); 
              end
           end
        end
        
        function r = rho(obj,t,type,varargin)
            %blocks of the second functions are embedded into
            %the blocks of the first function
            %one letter for the second, two letters for the first
            
            if ~isempty(varargin) %for p-p convention embed either in X_{p+1} or in X_p
                index = varargin{1};
            else
                index = obj.fun1.p; %in p-p-1, always embed into X_p
            end
           if strcmpi(type,'K')
              KK = obj.fun1.getInterval(index,'K');
              K = obj.fun2.getInterval(t,'K');
              if ~isempty(K) && ~isempty(KK)
                r = [KK(1)+(KK(2)-KK(1)+1)-(K(2)-K(1)+1),KK(2)];
              else
                r = cell(0);
              end
           elseif strcmpi(type,'L')
              KK = obj.fun1.getInterval(index,'L');
              K = obj.fun2.getInterval(t,'L');
              if ~isempty(K) && ~isempty(KK)
                r = [KK(1),KK(2)-(KK(2)-KK(1)+1)+(K(2)-K(1)+1)];
              else
                r = cell(0);
              end
           elseif strcmpi(type,'Phi')
              r = obj.rho(t,'K',index);
              Phi = obj.fun2.getInterval(t,'Phi');
              if ~isempty(r) && ~isempty(Phi)
                  r(2) = r(1)+Phi(2)-Phi(1); 
              else
                  r = cell(0);
              end
           elseif strcmpi(type,'Psi')
              r = obj.rho(t,'L',index);
              Psi = obj.fun2.getInterval(t,'Psi');
              if ~isempty(r) && ~isempty(Psi) 
                  r(1) = r(2)-(Psi(2)-Psi(1)); 
              else
                  r = cell(0);
              end
           elseif strcmpi(type,'K-Phi')
               Phi = obj.rho(t,'Phi',index);
               r = obj.rho(t,'K',index);
               if ~isempty(Phi) && ~isempty(r)
                   r(1) = Phi(2)+1;
                   if r(2) < r(1)
                      r = cell(0); 
                   end
               elseif isempty(Phi) && ~isempty(r)
                  return;
               else
                   r = cell(0);
                   return;
               end
           elseif strcmpi(type,'L-Psi')
               Psi = obj.rho(t,'Psi',index);
               r = obj.rho(t,'L',index);
               if ~isempty(Psi) && ~isempty(r)
                   r(2) = Psi(1)-1;
                   if r(2) < r(1)
                      r = cell(0); 
                   end
               elseif isempty(Psi) && ~isempty(r)
                  return;
               else
                   r = cell(0);
                   return;
               end
           end
        end
        
        function r = etll(obj)
            r = sym('0');
            p = obj.fun1.p;
           for t = 0:obj.fun2.nblocks+1
               b1 = obj.fun1.getBound(p,'b');
               b2 = obj.fun2.getBound(t,'b');
              if obj.compare(b1,'>',b2)
                  r = r+obj.b1(t)+obj.b2(t);
              elseif obj.compare(b1,'=',b2)
                  r = r+obj.b3(t); 
              end
           end
           r = r+obj.etal3;
        end
        
        function r = etrr(obj)
            %formula 4.51 on p.38 !without the diagonal part!
            
            hit = ones(1,3);
            r = sym('0');
            q = obj.fun1.q;
            ba1 = obj.fun1.getBound(obj.fun1.q,'ba');
           for t = 1:obj.fun2.nblocks
               ba2 = obj.fun2.getBound(t,'ba');
              if obj.compare(ba1,'>',ba2)
                  w1 = obj.bb1(t,'emb_index',q);
                  w2 = obj.bb2(t,'emb_index',q);
                  if w1 ~= sym('0') && hit(1) && obj.CALC_SHOW_MESSAGE
                      hit(1) = false;
                      fprintf('bb1 is hit; t = %d.\n',t);
                  end
                  if w2 ~= sym('0') && hit(2) && obj.CALC_SHOW_MESSAGE
                      hit(2) = false;
                      fprintf('bb2 is hit; t = %d.\n',t);
                  end
                  r = r+w1+w2;
              elseif obj.compare(ba1,'=',ba2)
                  w3 = obj.bb3(t,'emb_index',q);
                  r = r+w3; 
                  if w3 ~= sym('0') && hit(3) && obj.CALC_SHOW_MESSAGE
                      hit(3) = false;
                      fprintf('bb3 is hit; t = %d.\n',t);
                  end
              end
           end
           w4 = obj.etar4;
           if w4 ~= sym('0') && obj.CALC_SHOW_MESSAGE
              fprintf('Last sum is hit\n');
           end
           r = r+w4;
        end
        
        function r = xill(obj)
           
           r = sym('0');
           p = obj.fun1.p;
           hit = ones(1,5);
           %the first two sums
           b1 = obj.fun1.getBound(p,'b');
           
           if strcmpi(obj.CALC_CONVENTION,'p-p')
               Psip1 = obj.fun1.getInterval(p+1,'Psi');
               sigma_emb_index = p;
               bb1 = obj.fun1.getBound(p,'bb');
           else 
               Psip1 = obj.fun1.getInterval(p,'Psi');
               sigma_emb_index = p-1;
               bb1 = obj.fun1.getBound(p-1,'bb');
           end
           
           if strcmpi(obj.CALC_CONVENTION,'p-p')
              b11 = obj.fun1.getBound(p+1,'b');
              for t = 1:obj.fun2.nblocks
                  b2 = obj.fun2.getBound(t,'b');
                  if obj.compare(b11,'>=',b2)
                        w = obj.b2(t,'emb_index',p+1);
                      r = r+w;
                      if w ~= 0 && hit(1) && obj.CALC_SHOW_MESSAGE
                         fprintf('First sum is hit.\n');
                         hit(1) = 0;
                      end
                  end
               end
           end
            
           for t = 2:obj.fun2.nblocks+1
               b2 = obj.fun2.getBound(t,'b');
              if obj.compare(b1,'>=',b2)
                    w = obj.b2(t,'emb_index',p);
                  r = r+w;
                  if w ~= 0 && hit(1) && obj.CALC_SHOW_MESSAGE
                     fprintf('First sum is hit.\n');
                     hit(1) = 0;
                  end
              end
           end
           for t = 0:obj.fun2.nblocks
              bb2 = obj.fun2.getBound(t,'bb');
              condition = obj.compare(bb1,'<',bb2) || (isempty(Psip1) &&  obj.compare(bb1,'=',bb2));
              if condition%obj.compare(bb1,'<',bb2)
                  w = obj.bb4(t,'emb_index',sigma_emb_index);
                  r = r+w;
                  if w~=0 && hit(2) && obj.CALC_SHOW_MESSAGE
                      fprintf('Second sum is hit.\n');
                      hit(2) = 0;
                  end
              end
           end
           %the third term (which is the second row in the formula)
           if strcmpi(obj.CALC_CONVENTION,'p-p')
               bound = p+1;
           else
               bound = p;
           end
           for u=1:bound
               for t = 1:obj.fun2.nblocks
                   t1 = obj.move(obj.fun1.get(u,'grLL','L'),'J');
                   t2 = obj.Tinvc(obj.move(obj.fun2.get(t,'grLL','bL-Psi'),'bJ'));
                   w = trace(t1*t2);
                   r = r + w;
                  if w~=0 && hit(3) && obj.CALC_SHOW_MESSAGE
                     fprintf('Third sum is hit.\n');
                     hit(3) = 0;
                  end
               end
           end
           
           %the fourth term
           %this is redundant, but just to keep as close to the formula as
           %possible
           if strcmpi(obj.CALC_CONVENTION,'p-p')
               bound = p;
           else
               bound = p-1;
           end
           for u = 1:bound%p-1
               for t = 1:obj.fun2.nblocks
                    t1 = obj.move(obj.fun1.get(u,'grLL','bL-Psi'),'bJ');
                    t2 = obj.proj(obj.move(obj.fun2.get(t,'grLL','bL-Psi'),'bJ'), '2c');
                    w = trace(t1*t2);
                    r = r + w;
                    if w ~= 0 && hit(4) && obj.CALC_SHOW_MESSAGE
                         fprintf('Fourth sum is hit.\n');
                         hit(4) = 0;
                    end
               end
           end
           
           %the fifth term
           f1 = obj.s.get_f(obj.fun1.type,obj.fun1.i,obj.fun1.j);
           for t = 1:obj.fun2.nblocks
               t1 = obj.fun2.get(t,'grL','Psi+1','bK');
               t2 = obj.fun2.get(t,'L','bK','Psi+1');
               w = trace(t1*t2)*f1;
               if w ~= 0 && hit(5) && obj.CALC_SHOW_MESSAGE
                 fprintf('Fifth sum is hit.\n');
                 hit(5) = 0;
               end
              r = r + (obj.countb1(t)+obj.countb2(t))*w;
           end

        end
        
        function r = xirr(obj)
            %this realizes (4.69) in plethora
            
           r = sym('0');
           p = obj.fun1.p;
           hit = ones(1,6);
           
           %first term
           %it exists in Plethora's convention and doesn't exist in p-p
           if strcmpi(obj.CALC_CONVENTION,'p-p-1')
               ba1 = obj.fun1.getBound(p-1,'ba');
               for t = 1:obj.fun2.nblocks
                   ba2 = obj.fun2.getBound(t,'ba');
                  if obj.compare(ba2,'<=',ba1)
                      w = obj.bb2(t,'emb_index',p-1);
                      r = r + w;
                      if w ~= 0 && hit(1) && obj.CALC_SHOW_MESSAGE
                         fprintf('First sum is hit.\n');
                         hit(1) = 0;
                      end
                  end
               end
           end
           
           %second term
           ba1 = obj.fun1.getBound(p,'ba');
           for t = 1:obj.fun2.nblocks
               ba2 = obj.fun2.getBound(t,'ba');
              if obj.compare(ba2,'<=',ba1)
                  w = obj.bb2(t,'emb_index',p);
                  r = r + w;
                  if w ~= 0 && hit(2) && obj.CALC_SHOW_MESSAGE
                     fprintf('Second sum is hit.\n');
                     hit(2) = 0;
                  end
              end
           end

           %third term
           a1 = obj.fun1.getBound(p,'a');
           Phip = obj.fun1.getInterval(p,'Phi');
           isemptyPhip = isempty(Phip);
           for t = 1:obj.fun2.nblocks
               a2 = obj.fun2.getBound(t,'a');
               condition = obj.compare(a2,'>',a1);
               condition = condition || (obj.compare(a2,'==',a1)&& isemptyPhip);
              if condition%obj.compare(a2,'>',a1)
                  w = obj.b4(t);
                  r = r + w;
                  if w ~= 0 && hit(3) && obj.CALC_SHOW_MESSAGE
                     fprintf('Third sum is hit.\n');
                     hit(3) = 0;
                  end
              end
           end
           
           %fourth term
           %r = r + obj.xirr4();
           for u = 1:p
               for t = 1:obj.fun2.nblocks
                    t1 = obj.move(obj.fun1.get(u,'LgrL','bK'),'bI');
                    t2 = obj.Tr(obj.move(obj.fun2.get(t,'LgrL','K-Phi'), 'I'));
                    w = trace(t1*t2);
                    r = r + w;
                    if w ~= 0 && hit(4) && obj.CALC_SHOW_MESSAGE
                        fprintf('Fourth sum is hit.\n');
                        hit(4) = 0;
                    end
               end
           end
           
           %fifth term
           %r = r + obj.xirr5();
           for u = 1:p
               for t = 1:obj.fun2.nblocks
                    t1 = obj.move(obj.fun1.get(u,'LgrL','K-Phi'),'I');
                    t2 = obj.proj(obj.move(obj.fun2.get(t,'LgrL','K-Phi'), 'I'),'1r');
                    w = trace(t1*t2);
                    r = r + w;
                    if w ~= 0 && hit(5) && obj.CALC_SHOW_MESSAGE
                        fprintf('Fifth sum is hit.\n');
                        hit(5) = 0;
                    end
               end
           end
           
           %sixth term
           f1 = obj.s.get_f(obj.fun1.type,obj.fun1.i,obj.fun1.j);
           for t = 1:obj.fun2.nblocks
              t1 = obj.fun2.get(t,'L','Phi','L');
              t2 = obj.fun2.get(t,'grL','L','Phi');
              w = trace(t1*t2)*f1;
              if w ~= 0 && hit(6) && obj.CALC_SHOW_MESSAGE
                  fprintf('Sixth sum is hit.\n');
                  hit(6) = 0;
              end
              if w~=0
                 r = r + (obj.counta1(t)+obj.counta2(t))*w;
              end
           end
           
        end
        
        function r = xirr5(obj)
           p = obj.fun1.p;
           r = sym('0');
           for u = 1:p
               for t = 1:obj.fun2.nblocks
                    t1 = obj.move(obj.fun1.get(u,'LgrL','K-Phi'),'I');
                    t2 = obj.proj(obj.move(obj.fun2.get(t,'LgrL','K-Phi'), 'I'),'1r');
                    w = trace(t1*t2);
                    r = r + w;
               end
           end 
        end
        
        function r = xirr4(obj)
           r = sym('0');
           p = obj.fun1.p;
           for u = 1:p
               for t = 1:obj.fun2.nblocks
                    t1 = obj.move(obj.fun1.get(u,'LgrL','bK'),'bI');
                    t2 = obj.Tr(obj.move(obj.fun2.get(t,'LgrL','K-Phi'), 'I'));
                    w = trace(t1*t2);
                    r = r + w;
               end
           end 
        end
        
        function r = etal3(obj)
           %computes the third term (the last sum) in 4.36
           r = sym('0');
           p = obj.fun1.p;
           for t = 1:obj.fun2.nblocks+1
               b1 = obj.fun1.getBound(p,'b');
               b2 = obj.fun2.getBound(t,'b');
              if obj.compare(b1,'>',b2)
                  t1 = obj.fun1.geti('LgrL',obj.rho(t,'K'));
                  t2 = obj.fun2.get(t,'LgrL','K');
                  
                  t3 = obj.fun1.geti('grLL',obj.rho(t,'L'));
                  t4 = obj.fun2.get(t,'grLL','L');
                  
                  if obj.dimcheck2(t1,t2)&&obj.dimcheck2(t3,t4)
                       r = r + trace(t1*t2)-trace(t3*t4);
                  end
              end
           end
        end
        
        function r = etar4(obj)
           %computes the third term (the last sum) in 4.36
           r = sym('0');
           q = obj.fun1.q;
           ba1 = obj.fun1.getBound(q,'ba');
           for t = 1:obj.fun2.nblocks
               ba2 = obj.fun2.getBound(t,'ba');
              if obj.compare(ba1,'>',ba2)
                  t1 = obj.fun1.geti('grLL',obj.sigma(t,'bL',q));
                  t2 = obj.fun2.get(t,'grLL','bL');
                  
                  t3 = obj.fun1.geti('LgrL',obj.sigma(t,'bK',q));
                  t4 = obj.fun2.get(t,'LgrL','bK');
                  
                  if obj.dimcheck2(t1,t2) && obj.dimcheck2(t3,t4)
                       r = r + trace(t1*t2)-trace(t3*t4);
                  end
              end
           end
        end
        
        function r = b1(obj,t)
           t1 = obj.fun1.geti('LgrL',obj.rho(t,'Phi'));
           t2 = obj.fun2.get(t,'L','Phi','bL');
           t3 = obj.fun2.get(t,'grL','bL','Phi');
           if calcL.dimcheck(t1,t2,t3)
               r = -trace(t1*t2*t3);
           else
               r = 0;
           end
        end
        
        function r = b2(obj,t,varargin)
            replace = false;
            emb_index = obj.fun1.p;
            if ~isempty(varargin)
               for i = 1:length(varargin)
                  if strcmpi(varargin{i},'replace') 
                     if varargin{i+1}
                         replace = true;
                     end
                  elseif strcmpi(varargin{i},'emb_index')
                      emb_index = varargin{i+1};
                  end
               end
            end
            if replace
                interval = obj.fun1.getInterval(obj.fun1.p,'Psi');
            else
                interval = obj.rho(t,'Psi',emb_index);
            end
           t1 = obj.fun1.geti('grLL',interval);
           t2 = obj.fun2.get(t,'grL','Psi','bK-1');
           t3 = obj.fun2.get(t,'L','bK-1','Psi');
           if calcL.dimcheck(t1,t2,t3)
               r = trace(t1*t2*t3);
           else
               r = 0;
           end
        end
        
        function r = b3(obj,t)
           %version 1, with empty blocks as in plethora
%            p = obj.fun1.p;
%            t1 = obj.fun1.get(p,'grLL','Psi','L-Psi');
%            t2 = obj.fun2.get(t,'grL','L-Psi','K');
%            t3 = obj.fun2.get(t,'L','K','Psi');
%            if ~calcL.dimcheck(t1,t2,t3)
%                r = 0;
%                return;
%            end
%            r = trace(t1*t2*t3);
           
           %version 2, without empty blocks, as we want now:
           p = obj.fun1.p;
           rhoPsit = obj.rho(t,'Psi');
           Lp = obj.fun1.getInterval(p,'L');
           if ~isempty(Lp)
              diff = Lp;
              if ~isempty(rhoPsit)
                  diff(2) = rhoPsit(1)-1;
                  if diff(2) < diff(1)
                      diff = cell(0);
                  end
              end
           end
           t1 = obj.fun1.geti('grLL',rhoPsit,diff);
           t2 = obj.fun2.get(t,'grL','L-Psi','K');
           t3 = obj.fun2.get(t,'L','K','Psi');
           if ~calcL.dimcheck(t1,t2,t3)
               r = 0;
               return;
           end
           r = trace(t1*t2*t3);
        end
        
        function r = bb1(obj,t,varargin)
            emb_index = obj.fun1.q;
            if ~isempty(varargin)
               for i = 1:length(varargin)
                 if strcmpi(varargin{i},'emb_index')
                  emb_index = varargin{i+1};
                 end
               end
            end
           t1 = obj.fun1.geti('grLL',obj.sigma(t+1,'Psi',emb_index));
           t2 = obj.fun2.get(t+1,'grL','Psi','K');
           t3 = obj.fun2.get(t+1,'L','K','Psi');
           if ~calcL.dimcheck(t1,t2,t3)
               r = 0;
               return;
           end
           r = -trace(t1*t2*t3);
        end
        
        function r = bb2(obj,t,varargin)
            replace = false;
            num_replace = -1;
            emb_index = obj.fun1.p; %default in p-p convention
           if ~isempty(varargin)
              for i = 1:length(varargin)
                  if strcmpi(varargin{i},'replace')
                      replace = true;
                      num_replace = varargin{i+1};
                  elseif strcmpi(varargin{i},'emb_index')
                      %if we are in p-p-1 convention, might vary
                       emb_index = varargin{i+1};
                  end
              end
           end
           if ~replace
               interval = obj.sigma(t,'Phi',emb_index);
           elseif replace
               interval = obj.fun1.getInterval(num_replace,'Phi');
           end
           t1 = obj.fun1.geti('LgrL',interval);
           t2 = obj.fun2.get(t,'L','Phi','L');
           t3 = obj.fun2.get(t,'grL','L','Phi');
           if ~calcL.dimcheck(t1,t2,t3)
               r = 0;
               return;
           end
           r = trace(t1*t2*t3);
        end
        
        function r = bb3(obj,t,varargin)
            %OLD code that works with the convention on empty blocks
%            t1 = obj.fun1.get(obj.fun1.q,'LgrL','bK-Phi','Phi');
%            t2 = obj.fun2.get(t,'L','Phi','bL');
%            t3 = obj.fun2.get(t,'grL','bL','bK-Phi');
%            if ~calcL.dimcheck(t1,t2,t3)
%                r = 0;
%                return;
%            end
%            r = trace(t1*t2*t3);

            emb_index = obj.fun1.q;
            if ~isempty(varargin)
               for i = 1:length(varargin)
                 if strcmpi(varargin{i},'emb_index')
                   emb_index = varargin{i+1};
                 end
               end
            end
           %NEW code without the empty blocks
           sigmaPhit = obj.sigma(t,'Phi',emb_index);
           bKq = obj.fun1.getInterval(obj.fun1.q,'bK');
           %bKq = obj.fun1.getInterval(obj.fun1.p,'bK');
           diff = bKq;
           if ~isempty(bKq) && ~isempty(sigmaPhit)
               diff(2) = sigmaPhit(1)-1;
               if diff(2) < diff(1)
                   diff = cell(0);
               end
           end
           t1 = obj.fun1.geti('LgrL',diff,sigmaPhit);
           t2 = obj.fun2.get(t,'L','Phi','bL');
           t3 = obj.fun2.get(t,'grL','bL','bK-Phi');
           if ~calcL.dimcheck(t1,t2,t3)
               r = 0;
               return;
           end
           r = trace(t1*t2*t3);
        end
        
        function r = bb4(obj,t,varargin)
           emb_index = obj.fun1.p; %default in p-p convention
           if ~isempty(varargin)
              for i = 1:length(varargin)
                  if strcmpi(varargin{i},'emb_index')
                      %if we are in p-p-1 convention, might vary
                       emb_index = varargin{i+1};
                  end
              end
           end
           t1 = obj.fun1.geti('grLL',obj.sigma(t+1,'Psi',emb_index));
           t2 = obj.fun2.get(t,'grL','Psi+1','bK');
           t3 = obj.fun2.get(t,'L','bK','Psi+1');
           if ~calcL.dimcheck(t1,t2,t3)
               r = 0;
               return;
           end
           r = trace(t1*t2*t3); 
        end
        
        function r = b4(obj,t)
           t1 = obj.fun1.geti('LgrL',obj.rho(t,'Phi'));
           t2 = obj.fun2.get(t,'L','Phi','L');
           t3 = obj.fun2.get(t,'grL','L','Phi');
           if ~calcL.dimcheck(t1,t2,t3)
               r = 0;
               return;
           end
           r = trace(t1*t2*t3); 
        end
       
        function r = move(obj,A,type)
            %I realized this function only for square A's
            n = obj.s.n;
            r = sym(zeros(obj.s.n));
            [N,~] = size(A);
            if isempty(A)
                r = sym(zeros(n));
                return;
            end
           if strcmpi(type,'bJ')
               r(n-N+1:n, n-N+1:n) = A;
           elseif strcmpi(type,'J')
               r(1:N,1:N) = A;
           elseif strcmpi(type,'I')
               r(n-N+1:n, n-N+1:n) = A;
           elseif strcmpi(type,'bI')
               r(1:N, 1:N) = A;
           end
        end
        
        function s = countb1(obj,t)
            %counts the first term in the last sum of 4.63

            p = obj.fun1.p;
            s = 0;

            for u = 1:(p-1)
                b1 = obj.fun1.getBound(u,'b');
                b2 = obj.fun2.getBound(t+1,'b');%obj.fun2.getBound(t+1,'b');
                if obj.compare(b1,'>=',b2)
                   s = s + 1; 
                end
            end
        end
        
        function s = countb2(obj,t)
            %counts the second term in the last sum of 4.63
            p = obj.fun1.p;
            s = 0;
            if strcmpi(obj.CALC_CONVENTION,'p-p')
                bound = p-1;
            else
                bound = p-2;
            end
            for u = 1:bound%2:(p-1) 
                bb1 = obj.fun1.getBound(u,'bb');%obj.fun1.getBound(u-1,'bb');
                bb2 = obj.fun2.getBound(t,'bb');

                if obj.compare(bb1,'<',bb2)
                    s = s + 1; 
                end
            end
        end
        
        function s = counta1(obj,t,varargin)
           p = obj.fun1.p;
           s = 0;
           
           if strcmpi(obj.CALC_CONVENTION,'p-p')
               bound = p-1;
           elseif strcmpi(obj.CALC_CONVENTION, 'p-p-1')
               bound = p-2;
           end
           
           for u =1:bound%p-2%p-1%p-2
               ba1 = obj.fun1.getBound(u,'ba');
               ba2 = obj.fun2.getBound(t,'ba');

              if obj.compare(ba1,'>=',ba2)
                 s = s+1; 
              end
           end
        end
        
        function s = counta2(obj,t)
           p = obj.fun1.p;
           s = 0;
           for u = 1:p-1
               a1 = obj.fun1.getBound(u,'a');
               a2 = obj.fun2.getBound(t,'a');
               
              if obj.compare(a1,'<',a2)
                 s = s+1; 
              end
           end
        end
        
        function s = by1(obj)
           % the first sum in Lemma 4.17 page 60
           s = sym('0');
           ba1 = obj.fun1.getBound(obj.fun1.p-1,'ba');
           bb1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
           for t = 1:obj.fun2.nblocks+2
               ba2 = obj.fun2.getBound(t-1,'ba');
               bb2 = obj.fun2.getBound(t-1,'bb');
               if obj.compare(ba2,'<',ba1) && obj.compare(bb2,'>',bb1)
                  t1 = obj.fun1.geti('grLL', obj.sigma(t,'Psi',obj.fun1.p-1));
                  t2 = obj.fun2.get(t,'grLL','Psi');
                  if calcL.dimcheck2(t1,t2)
                     s = s + trace(t1*t2); 
                  end
               end
           end
        end
        
        function s = by2(obj)
           % the second sum in Lemma 4.17 page 60
           s = sym('0');
           ba1 = obj.fun1.getBound(obj.fun1.p-1,'ba');
           bb1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
           for t = 1:obj.fun2.nblocks+2
               ba2 = obj.fun2.getBound(t-1,'ba');
               bb2 = obj.fun2.getBound(t-1,'bb');
               if obj.compare(ba2,'~=',ba1) && obj.compare(bb2,'==',bb1)
                  t1 = obj.fun1.get(obj.fun1.p,'grLL', 'Psi');
                  t2 = obj.fun2.get(t,'grLL','Psi');
                  if calcL.dimcheck2(t1,t2)
                     s = s + trace(t1*t2); 
                  end
               end
           end
        end
        
        function s = by3(obj)
           s = sym('0');
           ba1 = obj.fun1.getBound(obj.fun1.p-1,'ba');
           bb1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
           for t = 1:obj.fun2.nblocks+2
               ba2 = obj.fun2.getBound(t-1,'ba');
               bb2 = obj.fun2.getBound(t-1,'bb');
               if obj.compare(ba2,'==',ba1) && obj.compare(bb2,'<',bb1)
                  t1 = obj.fun2.get(t-1,'L', 'bK','L');
                  t2 = obj.fun2.get(t-1,'grL','L','bK');
                  t3 = det(obj.fun1.fun);
                  if calcL.dimcheck2(t1,t2)
                     s = s + trace(t1*t2)*t3; 
                  end
               end
           end
        end
        
        function s = by4(obj)
           s = sym('0');
           ba1 = obj.fun1.getBound(obj.fun1.p-1,'ba');
           bb1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
           for t = 1:obj.fun2.nblocks+2
               ba2 = obj.fun2.getBound(t-1,'ba');
               bb2 = obj.fun2.getBound(t-1,'bb');
               if obj.compare(ba2,'==',ba1) && obj.compare(bb2,'>=',bb1)
                  t1 = obj.fun1.get(obj.fun1.p-1,'LgrL', 'bK');
                  t2 = obj.fun2.get(t-1,'LgrL','bK');
                  if calcL.dimcheck2(t1,t2)
                     s = s + trace(t1*t2);
                  end
               end
           end
        end
        
        function s = by5(obj)
           s = sym('0');
           ba1 = obj.fun1.getBound(obj.fun1.p-1,'ba');
           bb1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
           for t = 1:obj.fun2.nblocks+2
               ba2 = obj.fun2.getBound(t-1,'ba');
               bb2 = obj.fun2.getBound(t-1,'bb');
               if obj.compare(ba2,'==',ba1) && obj.compare(bb2,'>=',bb1)
                  t1 = obj.fun1.geti('grLL', obj.sigma(t-1,'bL-Psi',obj.fun1.p-1));
                  t2 = obj.fun2.get(t-1,'grLL','bL-Psi');
                  if calcL.dimcheck2(t1,t2)
                     s = s + trace(t1*t2);
                  end
               end
           end
        end
        
        function s = by6(obj)
           s = sym('0');
           ba1 = obj.fun1.getBound(obj.fun1.p-1,'ba');
           bb1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
           for t = 1:obj.fun2.nblocks+2
               ba2 = obj.fun2.getBound(t-1,'ba');
               bb2 = obj.fun2.getBound(t-1,'bb');
               if obj.compare(ba2,'==',ba1) && obj.compare(bb2,'==',bb1)...
                       && obj.compareExitY(obj.fun1.p-1,t-1)
                  t1 = obj.fun2.get(t-1,'L', 'Phi','L');
                  t2 = obj.fun2.get(t-1,'grL','L','Phi');
                  t3 = det(obj.fun1.fun);
                  if calcL.dimcheck2(t1,t2)
                     s = s + trace(t1*t2)*t3;
                  end
               end
           end
        end
        
        function s = by7(obj)
           s = sym('0');
           ba1 = obj.fun1.getBound(obj.fun1.p-1,'ba');
           bb1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
           for t = 1:obj.fun2.nblocks+2
               ba2 = obj.fun2.getBound(t-1,'ba');
               bb2 = obj.fun2.getBound(t-1,'bb');
               if obj.compare(ba2,'==',ba1) && obj.compare(bb2,'==',bb1)...
                       && obj.compareExitY(obj.fun1.p-1,t-1)
                  t1 = obj.fun1.get(obj.fun1.p-1,'grLL', 'bL');
                  t2 = obj.fun2.get(t-1,'grLL','bL');
                  if calcL.dimcheck2(t1,t2)
                     s = s + trace(t1*t2);
                  end
               end
           end 
        end
        
        function s = by8(obj)
           s = sym('0');
           ba1 = obj.fun1.getBound(obj.fun1.p-1,'ba');
           bb1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
           for t = 1:obj.fun2.nblocks+2
               ba2 = obj.fun2.getBound(t-1,'ba');
               bb2 = obj.fun2.getBound(t-1,'bb');
               if obj.compare(ba2,'==',ba1) && obj.compare(bb2,'==',bb1)...
                       && obj.compareExitY(obj.fun1.p-1,t-1)
                  t1 = obj.fun1.get(obj.fun1.p-1,'LgrL', 'bK');
                  t2 = obj.fun2.get(t-1,'LgrL','bK');
                  if calcL.dimcheck2(t1,t2)
                     s = s + trace(t1*t2);
                  end
               end
           end 
        end
        
        function s = by(obj)
            s = obj.by1 + obj.by2 + obj.by3 + obj.by4 + obj.by5 + ...
                obj.by6 + obj.by7 - obj.by8;
        end
        
        function s = bContribution(obj)
            %the contribution of all B-terms;
            %I want to compare it with Lemma 4.16 and Lemma 4.97 
            s = sym('0');
            p = obj.fun1.p;
            q = obj.fun1.q;
            
            b1 = obj.fun1.getBound(p,'b');
            %contributions from etaL:
           for t = 1:obj.fun2.nblocks
               b2 = obj.fun2.getBound(t,'b');
              if obj.compare(b1,'>',b2)
                  s = s-obj.b1(t)-obj.b2(t);
              elseif obj.compare(b1,'=',b2)
                  s = s-obj.b3(t); 
              end
           end
           
           %contributions from etaR:
           ba1 = obj.fun1.getBound(q,'ba');
           for t = 1:obj.fun2.nblocks
               ba2 = obj.fun2.getBound(t,'ba');
              if obj.compare(ba1,'>',ba2)
                  s = s-obj.bb1(t,'emb_index',q)-obj.bb2(t,'emb_index',q);
              elseif obj.compare(ba1,'=',ba2)
                  s = s-obj.bb3(t,'emb_index',q); 
              end
           end
           
           %contributions from xiL:
           b1 = obj.fun1.getBound(p,'b');
           bb1 = obj.fun1.getBound(p-1,'bb');
           Psip = obj.fun1.getInterval(p,'Psi');
           for t = 1:obj.fun2.nblocks+1
               b2 = obj.fun2.getBound(t,'b');
              if obj.compare(b1,'>=',b2)
                  w = obj.b2(t);
                  s = s+w;
              end
              bb2 = obj.fun2.getBound(t,'bb');
              condition = obj.compare(bb1,'<',bb2) || (obj.compare(bb1,'==',bb2)&&isempty(Psip));
              if condition%obj.compare(bb1,'<',bb2)
                  w = obj.bb4(t,'emb_index',p-1);
                  s = s+w;
              end
           end
           
           %contributions from xiR:
           %first term
           ba1 = obj.fun1.getBound(p-1,'ba');
           for t = 1:obj.fun2.nblocks
               ba2 = obj.fun2.getBound(t,'ba');
              if obj.compare(ba2,'<=',ba1)
                  w = obj.bb2(t,'emb_index',p-1);
                  s = s + w;
              end
           end
           
           %second term
           ba1 = obj.fun1.getBound(p,'ba');
           for t = 1:obj.fun2.nblocks
               ba2 = obj.fun2.getBound(t,'ba');
              if obj.compare(ba2,'<=',ba1)
                  w = obj.bb2(t,'emb_index',p);
                  s = s + w;
              end
           end

           %third term
           a1 = obj.fun1.getBound(p,'a');
           Phip = obj.fun1.getInterval(p,'Phi');
           for t = 1:obj.fun2.nblocks
               a2 = obj.fun2.getBound(t,'a');
               condition = obj.compare(a2,'>',a1) || (obj.compare(a2,'==',a1)&& isempty(Phip));
              if condition
                  w = obj.b4(t);
                  s = s + w;
              end
           end
           
        end
        
        function s = diffCompareB(obj)
           if strcmpi(obj.fun1.type,'h')
               s = obj.by - obj.bContribution;
           elseif strcmpi(obj.fun1.type,'g')
               %realize later
               s = -1;
           end
        end
        
        function b = compareExitY(obj,p,t)
           b = false;
           exit1 = obj.fun1.getBound(p,'exitY');
           exit2 = obj.fun2.getBound(t,'exitY');
           if exit1 == 0 || exit2 == 0
               return;
           else
              if exit2 < exit1 
                  b = true;
              end
           end
        end
        
        function b = compareExitX(obj,p,t)
           b = false;
           exit1 = obj.fun1.getBound(p,'exitX');
           exit2 = obj.fun2.getBound(t,'exitX');
           if exit1 == 0 || exit2 == 0
               return;
           else
              if exit2 < exit1 
                  b = true;
              end
           end
        end

    end
    methods (Access = private)
        function b = compare(obj,b1,sign,b2)
            if b1 == 0 || b2 == 0
                b = false;
                return;
            end
            if strcmpi(sign,'<')
                b = (b1 < b2);
            elseif strcmpi(sign,'<=')
                b = (b1 <= b2);
            elseif strcmpi(sign,'>')
                b = (b1 > b2);
            elseif strcmpi(sign,'=>') || strcmpi(sign,'>=')
                b = (b1 >= b2); 
            elseif strcmpi(sign,'==')||strcmpi(sign,'=')
                b = (b1 == b2); 
            elseif strcmpi(sign,'~=')
                b = (b1 ~= b2); 
            end
        end
    end
    methods
        %A list of methods for various tests
        %The enumeration follows Misha's Plethora
        %by default, the formulae are applied to obj.fun1
        %if need both, the second argument is the num (fun1 or fun2)
        
        %Section of formulae for xirr
        function A = formula472_1u(obj,u,varargin)
            %the teeh term in Formula 4.72(1)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            if num == 1
                A  = obj.proj(obj.move(obj.fun1.get(u,'LgrL','bK'),'bI'),'2r');
            elseif num == 2
                A  = obj.proj(obj.move(obj.fun2.get(u,'LgrL','bK'),'bI'),'2r');
            end
        end
        function A = formula472_1(obj,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            A = sym(zeros(obj.s.n));
            if num == 1
                 for u = 1:obj.fun1.nblocks
                    A  = A + obj.formula472_1u(u,1);
                 end
            elseif num == 2
             for u = 1:obj.fun2.nblocks
                A  = A + obj.formula472_1u(u,2);
             end 
            end
        end
        function A = formula472_2u(obj,u,varargin)
             num = 1;
             if ~isempty(varargin)
                num = varargin{1}; 
             end
             if num == 1
                A = obj.Tr(obj.move(obj.fun1.get(u,'LgrL','K-Phi'), 'I'));
             elseif num == 2
                A = obj.Tr(obj.move(obj.fun2.get(u,'LgrL','K-Phi'), 'I'));
             end
        end
        function A = formula472_2(obj,varargin)
             A = sym(zeros(obj.s.n));
             num = 1;
             if ~isempty(varargin)
                num = varargin{1}; 
             end
             if num == 1
                 for u = 1:obj.fun1.nblocks
                    A = A + obj.formula472_2u(u,1);
                 end
             elseif num == 2
                 for u = 1:obj.fun2.nblocks
                    A = A + obj.formula472_2u(u,2);
                 end
             end
        end
        function A = formula472(obj,varargin)
            %formula 472 is the expression for
            %$Pi_{\Gamma_2^r)(\eta_R^1)_{\geq 0}$
           num = 1;
           if ~isempty(varargin)
              num = varargin{1}; 
           end
           A  = obj.formula472_1(num) + obj.formula472_2(num);
        end
        
        function A = formula470_1t(obj,t,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            A = sym(zeros(obj.s.n));
            Z = A;
            if num == 1
                t1 = obj.fun1.get(t,'L','Phi','L');
                t2 = obj.fun1.get(t,'grL','L','Phi');
                where = obj.fun1.getRun(t,'Phi','X');
            elseif num == 2
                t1 = obj.fun2.get(t,'L','Phi','L');
                t2 = obj.fun2.get(t,'grL','L','Phi');
                where = obj.fun2.getRun(t,'Phi','X');
            end
            if ~isempty(where)
                Z(where(1):where(2),where(1):where(2)) = t1*t2;
                A = A + obj.Tr(Z);
            end
        end
        function A = formula470_1(obj,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            if num == 1
                numblocks = obj.fun1.nblocks;
            elseif num == 2
                numblocks = obj.fun2.nblocks; 
            end
            A = sym(zeros(obj.s.n));
            for t = 1:numblocks
                A = A + obj.formula470_1t(t,num);
            end 
        end
        function A = formula470_2t(obj,t,varargin)
            num = 1;
            if ~isempty(varargin)
                num = varargin{1};
            end
            if num == 1
                A = obj.Tr(obj.move(obj.fun1.get(t,'LgrL','K-Phi'),'I'));
            elseif num == 2
                A = obj.Tr(obj.move(obj.fun2.get(t,'LgrL','K-Phi'),'I'));
            end
        end
        function A = formula470_2(obj,varargin)
            A = sym(zeros(obj.s.n));
            num = 1;
            if ~isempty(varargin)
                num = varargin{1};
            end
            if num == 1
                for t = 1:obj.fun1.nblocks
                   A = A + obj.formula470_2t(t,num);
                end 
            elseif num == 2
                for t = 1:obj.fun2.nblocks
                   A = A + obj.formula470_2t(t,num);
                end 
            end
        end
        function A = formula470(obj,varargin)
            %that's a formula for 
            % $\gamma_r(X \nabla_X)
            num = 1;
            if ~isempty(varargin)
                num = varargin{1};
            end
            A = obj.formula470_1(num) + obj.formula470_2(num); 
        end
        
        %section of formulae for xill
        function A = formula464_1t(obj,t,varargin)
            %the teeth term in Formula 4.64(1) (I shift by 1 though)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            A = sym(zeros(obj.s.n));
            Z = A;
            if num == 1
                t1 = obj.fun1.get(t,'grL','Psi+1','bK');
                t2 = obj.fun1.get(t,'L','bK','Psi+1');
                where = obj.fun1.getRun(t+1,'Psi','Y');
            elseif num == 2
                t1 = obj.fun2.get(t,'grL','Psi+1','bK');
                t2 = obj.fun2.get(t,'L','bK','Psi+1');
                where = obj.fun2.getRun(t+1,'Psi','Y');
            end
            if ~isempty(where)
                Z(where(1):where(2),where(1):where(2)) = t1*t2;
                A = A + obj.Tinvc(Z);
            end
        end
        function A = formula464_2t(obj,t,varargin)
            num = 1;
            if ~isempty(varargin)
                num = varargin{1};
            end
            if num == 1
                A = obj.Tinvc(obj.move(obj.fun1.get(t,'grLL','bL-Psi'),'bJ'));
            elseif num == 2
                A = obj.Tinvc(obj.move(obj.fun2.get(t,'grLL','bL-Psi'),'bJ'));
            end
        end
        function A = formula464_1(obj,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            A = sym(zeros(obj.s.n));
            if num == 1
                numblocks = obj.fun1.nblocks; 
            elseif num == 2
                numblocks = obj.fun2.nblocks;
            end
            for t = 1:numblocks
                A = A + obj.formula464_1t(t,num); 
            end
        end
        function A = formula464_2(obj,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            A = sym(zeros(obj.s.n));
            if num == 1
                numblocks = obj.fun1.nblocks; 
            elseif num == 2
                numblocks = obj.fun2.nblocks;
            end
            for t = 1:numblocks
                A = A + obj.formula464_2t(t,num); 
            end
        end 
        function A = formula464(obj,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            A = obj.formula464_1(num) + obj.formula464_2(num);
        end
        
        function A = formula466_1u(obj,u,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            if num == 1
                A  = obj.proj(obj.move(obj.fun1.get(u,'grLL','L'),'J'),'1c');
            elseif num == 2
                A  = obj.proj(obj.move(obj.fun2.get(u,'grLL','L'),'J'),'1c');
            end
           %the code below makes the formula work in the old convention
           %what's the issue: we define L_u = Psi_u for the last u, so in the formula such
           %block doesn't shift to the location defined by the Y-run; the
           %code below shifts it to the right place
           if obj.s.empty_blocks_convention && u == obj.fun1.nblocks+1
               m = obj.fun1.get(u,'grLL','Psi');
               where = obj.fun1.getRun(u,'Psi','X');
               A = sym(zeros(obj.s.n));
               if ~isempty(where)
                   A(where(1):where(2), where(1):where(2)) = m;
               end
           end
        end
        function A = formula466_2u(obj,u,varargin)
            num = 1;
            if ~isempty(varargin)
                num = varargin{1};
            end
            if num == 1
                A = obj.Tinvc(obj.move(obj.fun1.get(u,'grLL','bL-Psi'),'bJ'));
            elseif num == 2
                A = obj.Tinvc(obj.move(obj.fun2.get(u,'grLL','bL-Psi'),'bJ'));
            end
        end
        function A = formula466_1(obj,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            if num == 1
                numblocks = obj.fun1.nblocks;
            elseif num == 2
                numblocks = obj.fun2.nblocks;
            end
            if obj.s.empty_blocks_convention
                numblocks = numblocks + 1;
            end
            A = sym(zeros(obj.s.n));
            for u = 1:numblocks
               A = A + obj.formula466_1u(u,num);
            end
        end
        function A = formula466_2(obj,varargin)
            num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
            if num == 1
                numblocks = obj.fun1.nblocks;
            elseif num == 2
                numblocks = obj.fun2.nblocks;
            end
            if obj.s.empty_blocks_convention
                numblocks = numblocks + 1;
            end
            A = sym(zeros(obj.s.n));
            for u = 1:numblocks
               A = A + obj.formula466_2u(u,num);
            end 
        end
        function A = formula466(obj,varargin)
           num = 1;
            if ~isempty(varargin)
               num = varargin{1}; 
            end
           A = obj.formula466_1(num) + obj.formula466_2(num);
        end
        function A = formula465(obj)
            A = sym(zeros(1,1));
            Z = A;
           for t = 1:obj.fun1.nblocks+1
              t1 = obj.fun1.get(t,'grL','Psi','K');
              t2 = obj.fun1.get(t,'L','K','Psi');
              if ~isempty(t1)&& ~isempty(t2)
                  t12 = t1*t2;
              else
                  t12 = Z;
              end
              t3 = obj.fun1.get(t,'grL','Psi','bK-1');
              t4 = obj.fun1.get(t,'L','bK-1','Psi');
              if ~isempty(t3) && ~isempty(t4)
                  t34 = t3*t4;
              else
                 t34 = Z;
              end
              t5 = obj.fun1.get(t,'grLL','Psi');
              if isempty(t5)
                  t5 = Z;
              end
              if isa(t12+t34-t5,'sym')
                A = A + norm(expand(t12+t34-t5));
              end
           end
        end
        
        function A = test(obj,varargin)
%            r = sym('0');
%            for t = 1:obj.fun1.nblocks
%                r = r + obj.Tr(obj.move(obj.fun1.get(t,'LgrL','K-Phi'), 'I'));
%            end

            %formula 4.64
%             A = sym(zeros(obj.s.n));
%             Z = A;
%             for t = 2:obj.fun1.nblocks+1
%                t1 = obj.fun1.get(t-1,'grL','Psi+1','bK');
%                t2 = obj.fun1.get(t-1,'L','bK','Psi+1');
%                where = obj.fun1.getRun(t,'Psi','Y');
%                if ~isempty(where)
%                  Z(where(1):where(2),where(1):where(2)) = t1*t2;
%                  A = A + obj.Tinvc(Z);
%                end
%                A = A + obj.Tinvc(obj.move(obj.fun1.get(t-1,'grLL','bL-Psi'),'bJ'));
%             end

             %formula 4.70
%             A = sym(zeros(obj.s.n));
%             Z = A;
%             for t = 1:obj.fun1.nblocks
%                t1 = obj.fun1.get(t,'L','Phi','L');
%                t2 = obj.fun1.get(t,'grL','L','Phi');
%                where = obj.fun1.getRun(t,'Phi','X');
%                if ~isempty(where)
%                  Z(where(1):where(2),where(1):where(2)) = t1*t2;
%                  A = A + obj.Tr(Z);
%                end
%                A = A + obj.Tr(obj.move(obj.fun1.get(t,'LgrL','K-Phi'),'I'));
%             end

             %formula 4.66 for etaL
%              A = sym(zeros(obj.s.n));
%              for t = 1:obj.fun1.nblocks+1
%                 A = A + (obj.move(obj.fun1.get(t,'grLL','L'), 'J'));
%                 %if t > 1
%                     A  = A + obj.Tinvc(obj.move(obj.fun1.get(t,'grLL','bL-Psi'),'bJ'));
%                 %end
%              end
           
             
             if ~isempty(varargin)
                 for i = 1:length(varargin)
                     if strcmpi(varargin{i},'Log')
                        A = simplifyFraction(A/det(obj.fun1.fun),'Expand',true);
                     end
                     if strcmpi(varargin{i},'Proj')||strcmpi(varargin{i},'Diag')
                        A = obj.proj_d(A);
                     end
                 end
             end
             
        end
        
        function A = test4702(obj)
           %this tests the following sentence on page 45:
           %Therefore, the contribution of the second term in (4.70)
           %to the final result equals
           A = trace(obj.formula472(1)*obj.formula470_2(2));
           B = obj.xirr4+obj.xirr5;
           A = A-B;
        end 
        
        function A = testb4(obj)
           %according to the text on page 46,
           %b4 = <4.72(2),4.70(1)>, under the condition
           %\alpha_p^1 < \alpha_t^2
           %for \alpha_p^1 >= \alpha_t^2, must vanish
           A = trace(obj.formula472_2(1)*obj.formula470_1(2));
           for t = 1:obj.fun2.nblocks
               if obj.compare(obj.fun1.getBound(obj.fun1.p,'a'),'<',obj.fun2.getBound(t,'a'))
                   A = A - obj.b4(t);
               end
           end
        end
        
        function A = testbb2tu(obj,t,u)
           %tests the second term in 4.69 
           %A = trace(obj.formula472_1(1)*obj.formula470_1(2));
           A = 0;
           if obj.compare(obj.fun1.getBound(u,'ba'),'>',obj.fun2.getBound(t,'ba'))
               A = trace(obj.formula472_1u(u,1)*obj.formula470_1t(t,2)) - obj.bb2(t);
           elseif obj.compare(obj.fun1.getBound(u,'ba'),'=',obj.fun2.getBound(t,'ba'))
               A = trace(obj.formula472_1u(u,1)*obj.formula470_1t(t,2)) - obj.bb2(t,'replace',true);
           end
        end
        
        function [actual,current,mine,diff1,diff2] = testAlek(obj)
            %Alek asked me whether in the old convention there's 
            %a counterexample that shows we need an extra term in (4.63)
            %I use this function to find a counterexample
            actual = trace(obj.formula466_2(1)*obj.formula464_1(2));
            current = 0;
            mine = 0;
            psip = isempty(obj.fun1.getInterval(obj.fun1.p,'Psi'));
            bp1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
            for t = 1:obj.fun2.nblocks+1
                bt = obj.fun2.getBound(t,'bb');
                if obj.compare(bp1,'<',bt)
                    term = obj.bb4(t,'emb_index',obj.fun1.p-1);
                    current = current + term;
                    mine = mine + term;
                end
                if psip && obj.compare(bp1,'==',bt)
                   term = obj.bb4(t,'emb_index',obj.fun1.p-1);
                   mine = mine + term;
                end
            end
            diff1 = actual - current;
            diff2 = actual - mine;
        end
        
        function [actual,current,mine,diff1,diff2] = testAlek2(obj)
            %this test is to show Alek that we need 2 sigmas
            %upd: a counterexample was found in:
            %testplet(1,4,[1 4], [4 2], 5); <g21,g51>
            actual = trace(obj.formula466_2(1)*obj.formula464_1(2));
            current = 0;
            mine = 0;
            bp1 = obj.fun1.getBound(obj.fun1.p-1,'bb');
            for t = 1:obj.fun2.nblocks+1
                bt = obj.fun2.getBound(t,'bb');
                if obj.compare(bp1,'<',bt)
                    termq = obj.bb4(t,'emb_index',obj.fun1.q);
                    current = current + termq;
                    termp = obj.bb4(t,'emb_index',obj.fun1.p-1);
                    mine = mine + termp;
                end
            end
            diff1 = actual - current;
            diff2 = actual - mine;
        end
    end
    
    methods
       %a list of methods to work with the contribution of all B-terms
       %in case X is the leading block
       %the methods compute the difference between formulae
       %I used the other methods to check intermediate formulae
       function [r1,r2,diff] = bx1(obj)
           ap = obj.fun1.getBound(obj.fun1.p,'a');
           bp = obj.fun1.getBound(obj.fun1.p,'b');
           r1 = sym('0');
           r2 = sym('0');
           for t = 1:obj.fun2.nblocks
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'>',ap) && obj.compare(bt,'<',bp)
                   t1 = obj.fun1.geti('LgrL',obj.rho(t,'Phi'));
                   t2 = obj.fun2.get(t,'LgrL','Phi');
                   if calcL.dimcheck2(t1,t2)
                       t12 = trace(t1*t2);
                   else
                       t12 = sym('0');
                   end
                   r1 = r1 + obj.b4(t)-obj.b1(t);
                   r2 = r2 + t12;
               end
           end
           diff = r1-r2;
       end
       function [r1,r2,diff] = bx2(obj)
           ap = obj.fun1.getBound(obj.fun1.p,'a');
           bp = obj.fun1.getBound(obj.fun1.p,'b');
           r1 = sym('0');
           r2 = sym('0');
           phip = isempty(obj.fun1.getInterval(obj.fun1.p,'Phi'));
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'==',ap) && obj.compare(bt,'<',bp)
                   %t1 = obj.fun1.get(obj.fun1.p,'LgrL','Phi'); %old conv
                   t1 = obj.fun1.geti('LgrL',obj.rho(t,'Phi')); %new conv
                   t2 = obj.fun2.get(t,'LgrL','Phi');
                   if calcL.dimcheck2(t1,t2)
                       t12 = trace(t1*t2);
                   else
                       t12 = sym('0');
                   end
                   r1 = r1 + obj.bb2(t)-obj.b1(t);
                   if phip
                       r1 = r1 + obj.b4(t);
                   end
                   r2 = r2 + t12;
               end
           end
           diff = r1-r2;
       end
       function [r1,r2,diff] = bx3(obj)
           ap = obj.fun1.getBound(obj.fun1.p,'a');
           bp = obj.fun1.getBound(obj.fun1.p,'b');
           r1 = sym('0');
           r2 = sym('0');
           phip = isempty(obj.fun1.getInterval(obj.fun1.p,'Phi'));
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'==',ap) && obj.compare(bt,'>',bp)
                   %t1 = obj.fun1.get(obj.fun1.p,'LgrL','Phi'); %old conv
                   t1 = obj.fun1.geti('LgrL',obj.rho(t,'Phi')); %new conv
                   t2 = obj.fun2.get(t,'LgrL','Phi');
                   if calcL.dimcheck2(t1,t2)
                       t12 = trace(t1*t2);
                   else
                       t12 = sym('0');
                   end
                   r1 = r1 + obj.bb2(t);
                   if phip
                       r1 = r1 + obj.b4(t);
                   end
                   r2 = r2 + t12;
               end
           end
           diff = r1-r2;
       end
       function [r1,r2,diff] = bx4(obj)
           ap = obj.fun1.getBound(obj.fun1.p,'a');
           bp = obj.fun1.getBound(obj.fun1.p,'b');
           r1 = sym('0');
           r2 = sym('0');
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'<',ap) && obj.compare(bt,'==',bp)
                   %t1 = obj.fun1.get(obj.fun1.p,'LgrL','Phi'); %old conv
                   t1 = obj.fun2.get(t,'L','bK-1','L');
                   t2 = obj.fun2.get(t,'grL','L','bK-1');
                   if calcL.dimcheck2(t1,t2)
                      t12 = trace(t1*t2)*det(obj.fun1.fun); 
                   else
                       t12 = sym('0');
                   end
                   r1 = r1 + obj.b2(t)-obj.b3(t);
                   r2 = r2 + t12;
               end
           end
           diff = r1-r2;
       end
       function [r1,r2,diff] = bx5(obj)
           ap = obj.fun1.getBound(obj.fun1.p,'a');
           bp = obj.fun1.getBound(obj.fun1.p,'b');
           r1 = sym('0');
           r2 = sym('0');
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'>',ap) && obj.compare(bt,'==',bp)
                   t1 = obj.fun1.get(obj.fun1.p,'grLL','L');
                   t2 = obj.fun2.get(t,'grLL','L');
                   if calcL.dimcheck2(t1,t2)
                      t12 = trace(t1*t2); 
                   else
                       t12 = sym('0');
                   end
                   t3 = obj.fun1.geti('LgrL',obj.rho(t,'K-Phi'));
                   t4 = obj.fun2.get(t,'LgrL','K-Phi');
                   if calcL.dimcheck2(t3,t4)
                      t34 = trace(t3*t4); 
                   else
                       t34 = sym('0');
                   end
                   r1 = r1 + obj.b2(t)-obj.b3(t)+obj.b4(t);
                   r2 = r2 + t12-t34;
               end
           end
           diff = r1-r2;
       end
       function [r1,r2,diff] = bx6(obj)
           ap = obj.fun1.getBound(obj.fun1.p,'a');
           bp = obj.fun1.getBound(obj.fun1.p,'b');
           r1 = sym('0');
           r2 = sym('0');
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'==',ap) && obj.compare(bt,'==',bp)
                   comparex = obj.compareExitX(obj.fun1.p,t);
                   t1 = obj.fun1.get(obj.fun1.p,'grLL','L');
                   t2 = obj.fun2.get(t,'grLL','L');
                   if calcL.dimcheck2(t1,t2)
                      t12 = trace(t1*t2); 
                   else
                       t12 = sym('0');
                   end
                   t3 = t1;
                   t4 = t2;
                   if calcL.dimcheck2(t3,t4) && comparex
                       t34 = -trace(t3*t4); 
                   else
                       t34 = sym('0');
                   end
                   t5 = obj.fun1.geti('LgrL',obj.rho(t,'K-Phi'));
                   t6 = obj.fun2.get(t,'LgrL','K-Phi');
                   if calcL.dimcheck2(t5,t6)
                      t56 = -trace(t5*t6); 
                   else
                       t56 = sym('0');
                   end
                   t7 = obj.fun2.get(t,'L','bK-1','Psi');
                   t8 = obj.fun2.get(t,'grL','Psi','bK-1');
                   if calcL.dimcheck2(t7,t8) && comparex
                      t78 = trace(t7*t8)*det(obj.fun1.fun); 
                   else
                       t78 = sym('0');
                   end
                   t9 = obj.fun1.get(obj.fun1.p,'LgrL','K');
                   t10 = obj.fun2.get(t,'LgrL','K');
                   if calcL.dimcheck2(t9,t10) && comparex
                      t910 = trace(t9*t10); 
                   else
                      t910 = sym('0'); 
                   end
                   r1 = r1 + obj.b2(t)-obj.b3(t)+obj.b4(t);
                   r2 = r2 + t12 + t34 + t56 + t78 + t910;
                   
                   %test on extra
%                    [thp,~,~,~] = obj.thetas(t);
%                    t11 = obj.fun1.geti('lgrl',obj.fun1.getInterval(obj.fun1.p,'K'),thp);
%                    t12 = obj.fun1.geti('L',thp,obj.rho(t,'L-Psi'));
%                    t13 = obj.fun2.get(t,'grL','L-Psi','K');
%                    if calcL.dimcheck(t11,t12,t13)
%                       extra = -trace(t11*t12*t13); 
%                    else
%                        extra = 0;
%                    end
%                    r2 = r2 + extra;
                   
                   %r2 = r2 + t12 + t56;
                   if obj.CALC_SHOW_MESSAGE
                      if t12 ~= 0
                          fprintf('Hit <grLL,grLL>;\n');
                      end
                      if t34 ~=0
                          fprintf('Hit <grLL,grLL>^a;\n');
                      end
                      if t56 ~=0
                          fprintf('Hit <LgrL,LgrL>;\n');
                      end
                      if t78 ~=0
                          fprintf('Hit <L,grL>^a fun1;\n');
                      end
                      if t910 ~=0
                         fprintf('Hit <LgrL,LgrL>^a>;\n'); 
                      end
                   end
               end
           end
           diff = r1-r2;
       end
       function [r1,r2,diff] = bx7(obj)
           bp = obj.fun1.getBound(obj.fun1.p-1,'bb');
           r1 = sym('0');
           r2 = sym('0');
           for t = 0:obj.fun2.nblocks
               bt = obj.fun2.getBound(t,'bb');
               if obj.compare(bt,'>',bp)
                   %t1 = obj.fun1.get(obj.fun1.p,'LgrL','Phi'); %old conv
                   t1 = obj.fun2.get(t,'L','bK','Psi+1');
                   t2 = obj.fun2.get(t,'grL','Psi+1','bK');
                   if calcL.dimcheck2(t1,t2)
                      t12 = trace(t1*t2)*det(obj.fun1.fun); 
                   else
                       t12 = sym('0');
                   end
                   r1 = r1 + obj.bb4(t,'emb_index',obj.fun1.p-1);
                   r2 = r2 + t12;
               end
           end
           diff = r1-r2;
       end
       function [r1,r2,diff] = bx8(obj)
           ap = obj.fun1.getBound(obj.fun1.p-1,'ba');
           r1 = sym('0');
           r2 = sym('0');
           for t = 0:obj.fun2.nblocks
               at = obj.fun2.getBound(t,'ba');
               if obj.compare(at,'<=',ap)
                   %t1 = obj.fun1.get(obj.fun1.p,'LgrL','Phi'); %old conv
                   t1 = obj.fun2.get(t,'L','Phi','L');
                   t2 = obj.fun2.get(t,'grL','L','Phi');
                   if calcL.dimcheck2(t1,t2)
                      t12 = trace(t1*t2)*det(obj.fun1.fun); 
                   else
                       t12 = sym('0');
                   end
                   r1 = r1 + obj.bb2(t,'emb_index',obj.fun1.p-1);
                   r2 = r2 + t12;
               end
           end
           diff = r1-r2;
       end
       
       function [r1,r2,diff] = quicktest(obj)
           r1 = 0;
           r2 = 0;
          for t = 1:obj.fun1.nblocks
             t1 = obj.fun1.get(t,'grL','L','K');
             t2 = obj.fun1.get(t,'L','K','L-Psi');
             tt = obj.fun1.get(t,'grLL','L','L-Psi');
             r1 = r1+t1*t2;
             if ~isempty(tt)
                 r2 = r2 + tt;
             end
          end
          diff = r1-r2;
       end
       
       function [thp,tht,xip,xit] = thetas(obj,t)
           %this code determines ThetaP and ThetaT 
           %only need to catch a silly mistake in the proof
            p = obj.fun1.p;
            minpt = min(p,t);
            thp = cell(0);
            xip = cell(0);
            tht = cell(0);
            xit = cell(0);
            
            for k = 1:minpt                
                ba1 = obj.fun1.getBound(p-k,'ba');
                bb1 = obj.fun1.getBound(p-k,'bb');
                ba2 = obj.fun2.getBound(t-k,'ba');
                bb2 = obj.fun2.getBound(t-k,'bb');
                
                if obj.compare(ba1,'==',ba2) && obj.compare(bb1,'==',bb2)
                   a1 = obj.fun1.getBound(p-k,'a');
                   b1 = obj.fun1.getBound(p-k,'b');
                   a2 = obj.fun2.getBound(t-k,'a');
                   b2 = obj.fun2.getBound(t-k,'b');
                   if ~( obj.compare(a1,'==',a2) && obj.compare(b1,'==',b2)) && a1*a2*b1*b2 ~= 0
                      %covers ia and ib, basically
                      [thp,tht,xip,xit] = obj.thetasAssign(t,k); 
                      break;
                   elseif a1*a2*b1*b2 == 0
                      [thp,tht,xip,xit] = obj.thetasAssign(t,k); 
                      break;
                   end
                elseif ba1*bb1*ba2*bb2 ~= 0
                    %found non-coinciding non-trivial Y-blocks -> assign
                    [thp,tht,xip,xit] = obj.thetasAssign(t,k);
                    break;
                elseif ba1*bb1*ba2*bb2 == 0
                    [thp,tht,xip,xit] = obj.thetasAssign(t,k);
                end
            end           
       end

       function [thp,tht,xip,xit] = thetasAssign(obj,t,k)
          %an intermediate method for obj.thetas 
            thp = cell(0);
            tht = cell(0);
            xip = cell(0);
            xit = cell(0);
            
            p = obj.fun1.p;
            Kp = obj.fun1.getInterval(p,'K');
            Kt = obj.fun2.getInterval(t,'K');
            bLp = obj.fun1.getInterval(p-1,'bL');
            bLt = obj.fun2.getInterval(t-1,'bL');
            bKpm = obj.fun1.getInterval(p-k,'bK');
            bKtm = obj.fun2.getInterval(t-k,'bK');
            bLpm = obj.fun1.getInterval(p-k,'bL');
            bLtm = obj.fun2.getInterval(t-k,'bL');
            
            if ~isempty(Kp) && ~isempty(bKpm)
                thp = zeros(1,2);
                thp(1) = Kp(2) + 1;
                thp(2) = bKpm(2);
            end
            
            if ~isempty(Kt) && ~isempty(bKtm)
                tht = zeros(1,2);
                tht(1) = Kt(2) + 1;
                tht(2) = bKtm(2);
            end
            
            if ~isempty(bLp) && ~isempty(bLpm)
                xip = zeros(1,2);
                xip(1) = bLp(1);
                xip(2) = bLpm(2)-1;
                if xip(1) > xip(2)
                    xip = cell(0);
                end
            end
            
            if ~isempty(bLt) && ~isempty(bLtm)
                xit = zeros(1,2);
                xit(1) = bLt(1);
                xit(2) = bLtm(2)-1;
                if xit(1) > xit(2)
                    xit = cell(0);
                end
            end

       end
       
       function [r1,r2,diff] = bx6start(obj)
           %just to make sure
           %this gives the formula that branches into 8 cases with exit points in bx6
           %it works, sure
           p = obj.fun1.p;
           
           ap = obj.fun1.getBound(p,'a');
           bp = obj.fun1.getBound(p,'b');
           r1 = sym('0');
           r2 = sym('0');
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'==',ap) && obj.compare(bt,'==',bp)
                    t1 = obj.fun1.geti('grLL', obj.rho(t,'Psi'));
                    t2 = obj.fun2.get(t,'grLL','Psi');
                    if calcL.dimcheck2(t1,t2)
                        t12 = trace(t1*t2);
                    else
                        t12 = 0;
                    end
                    t3 = obj.fun1.geti('grLL',obj.rho(t,'Psi'), obj.fun1.getInterval(p,'L'));
                    t4 = obj.fun2.get(t,'grL','L','K');
                    t5 = obj.fun2.get(t,'L','K','Psi');
                    if calcL.dimcheck(t3,t4,t5)
                       t345 = -trace(t3*t4*t5); 
                    else
                        t345 = 0;
                    end
                    t6 = obj.fun1.geti('LgrL',obj.rho(t,'Phi'));
                    t7 = obj.fun2.get(t,'L','Phi','L');
                    t8 = obj.fun2.get(t,'grL','L','Phi');
                    if calcL.dimcheck(t6,t7,t8)
                        t678 = trace(t6*t7*t8);
                    else
                        t678 = 0;
                    end
                    r1 = r1 + obj.b2(t)-obj.b3(t)+obj.b4(t);
                    r2 = r2 + t12+t345+t678;
               end
           end
           diff = r1-r2;
       end
       
       function [r1,r2,diff] = bx6start2(obj)
           %the formula before we run into the cases
           p = obj.fun1.p;
           ap = obj.fun1.getBound(p,'a');
           bp = obj.fun1.getBound(p,'b');
           r1 = sym('0');
           r2 = sym('0');
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'==',ap) && obj.compare(bt,'==',bp)
                    t1 = obj.fun1.get(p,'grLL', 'L');
                    t2 = obj.fun2.get(t,'grLL','L');
                    if calcL.dimcheck2(t1,t2)
                        t12 = trace(t1*t2);
                    else
                        t12 = 0;
                    end
                    t3 = obj.fun1.geti('LgrL',obj.rho(t,'K-Phi'));
                    t4 = obj.fun2.get(t,'LgrL','K-Phi');
                    if calcL.dimcheck2(t3,t4)
                       t345 = -trace(t3*t4); 
                    else
                        t345 = 0;
                    end
                    [thp,~,~,~] = obj.thetas(t);
                    t6 = obj.fun1.geti('LgrL',obj.fun1.getInterval(p,'K'), thp);
                    %t7 = obj.fun1.geti('L',thp,obj.rho(t,'Psi'));
                    %t8 = obj.fun2.get(t,'grL','Psi','K');
                    t7 = obj.fun1.geti('L',thp,obj.fun1.getInterval(p,'L'));
                    t8 = obj.fun2.get(t,'grL','L','K');
                    if calcL.dimcheck(t6,t7,t8)
                        t678 = -trace(t6*t7*t8);
                    else
                        t678 = 0;
                    end
                    r1 = r1 + obj.b2(t)-obj.b3(t)+obj.b4(t);
                    r2 = r2 + t12+t345+t678;
                    
%                    t11 = obj.fun1.geti('LgrL',obj.fun1.getInterval(obj.fun1.p,'K'),thp);
%                    t12 = obj.fun1.geti('L',thp,obj.rho(t,'L-Psi'));
%                    t13 = obj.fun2.get(t,'grL','L-Psi','K');
%                    if calcL.dimcheck(t11,t12,t13)
%                       extra = -trace(t11*t12*t13); 
%                    else
%                        extra = 0;
%                    end
%                    r2 = r2 + extra;
               end
           end
           diff = r1-r2;
       end
       
       function [r1,r2,diff] = bx6start3(obj)
           %the formula before we run into the cases
           p = obj.fun1.p;
           ap = obj.fun1.getBound(p,'a');
           bp = obj.fun1.getBound(p,'b');
           r1 = sym('0');
           r2 = sym('0');
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'==',ap) && obj.compare(bt,'==',bp)
                    t1 = obj.fun1.get(p,'grLL', 'L');
                    t2 = obj.fun2.get(t,'grLL','L');
                    if calcL.dimcheck2(t1,t2)
                        t12 = trace(t1*t2);
                    else
                        t12 = 0;
                    end
                    t3 = obj.fun1.geti('LgrL',obj.rho(t,'Phi'));
                    t4 = obj.fun2.get(t,'L','Phi','L');
                    t5 = obj.fun2.get(t,'grL','L','Phi');
                    if calcL.dimcheck(t3,t4,t5)
                       t345 = trace(t3*t4*t5); 
                    else
                        t345 = 0;
                    end
                    [thp,~,~,~] = obj.thetas(t);
                    Kp = obj.fun1.getInterval(p,'K');
                    Lp = obj.fun1.getInterval(p,'L');
                    interval = unis(thp,Kp);
                    t6 = obj.fun1.geti('LgrL',Kp, interval);
                    t7 = obj.fun1.geti('L',interval,Lp);
                    t8 = obj.fun2.get(t,'grL','L','K');
                    if calcL.dimcheck(t6,t7,t8)
                        t678 = -trace(t6*t7*t8);
                    else
                        t678 = 0;
                    end
                    r1 = r1 + obj.b2(t)-obj.b3(t)+obj.b4(t);
                    r2 = r2 + t12+t345+t678;
               end
           end
           diff = r1-r2;
       end
       
       function [r1,r2,diff] = bx6start4(obj)
           %the formula before we run into the cases
           p = obj.fun1.p;
           ap = obj.fun1.getBound(p,'a');
           bp = obj.fun1.getBound(p,'b');
           r1 = sym('0');
           r2 = sym('0');
           for t = 1:obj.fun2.nblocks+1
               at = obj.fun2.getBound(t,'a');
               bt = obj.fun2.getBound(t,'b');
               if obj.compare(at,'==',ap) && obj.compare(bt,'==',bp)
                    t1 = obj.fun1.get(p,'grLL', 'L');
                    t2 = obj.fun2.get(t,'grLL','L');
                    if calcL.dimcheck2(t1,t2)
                        t12 = trace(t1*t2);
                    else
                        t12 = 0;
                    end
                    t3 = obj.fun1.geti('LgrL',obj.fun1.getInterval(p,'K'), obj.rho(t,'K-Phi'));
                    t4 = obj.fun1.geti('L',obj.rho(t,'K-Phi'),obj.fun1.getInterval(p,'L'));
                    t5 = obj.fun2.get(t,'grL','L','K');
                    if calcL.dimcheck(t3,t4,t5)
                       t345 = -trace(t3*t4*t5); 
                    else
                        t345 = 0;
                    end
                    [thp,~,~,~] = obj.thetas(t);
                    t6 = obj.fun1.geti('LgrL',obj.fun1.getInterval(p,'K'), thp);
                    t7 = obj.fun1.geti('L',thp,obj.fun1.getInterval(p,'L'));
                    t8 = obj.fun2.get(t,'grL','L','K');
                    if calcL.dimcheck(t6,t7,t8)
                        t678 = -trace(t6*t7*t8);
                    else
                        t678 = 0;
                    end
                    r1 = r1 + obj.b2(t)-obj.b3(t)+obj.b4(t);
                    r2 = r2 + t12+t345+t678;
               end
           end
           diff = r1-r2;
       end
       
    end
    
    methods (Static)
        function r = dimcheck2(t1,t2)
           if isempty(t1)||isempty(t2)
               r = false;
               return;
           end
           [n1,m1] = size(t1);
           [n2,m2] = size(t2);
           if n1~= m2 || m1 ~= n2
               r = false;
           else
               r = true;
           end
        end
        function r = dimcheck(t1,t2,t3)
           %if the dimensions agree for matrix multiplication, r = true;
           %else it's false
           if isempty(t1)||isempty(t2)||isempty(t3)
               r = false;
               return;
           end
           [n1,m1] = size(t1);
           [n2,m2] = size(t2);
           [n3,m3] = size(t3);
           if m1 ~= n2 || m2 ~= n3 || n1 ~= m3
               r = false;
           else
               r = true;
           end
        end
    end
    
end

