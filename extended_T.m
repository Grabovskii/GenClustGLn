classdef extended_T < BD_data
    %this class realizes methods that allow to work with T (or gamma circ in
    %Misha's notation) on the level of Lie algebras
    
    %this version of T simply permutes gln-blocks.
    
    properties
        G1_pattern_r, G2_pattern_r %the way the projection from gln to glG1 works
        G1_pattern_c, G2_pattern_c %the way the projection from gln to glG1 works
        XYruns_r %an array of all non-trivial runs
        XYruns_c
    end
    
    methods
        function obj = extended_T(G1_r,G2_r,G1_c,G2_c,n)
            if nargin == 3
                n = G1_c;
                G1_c = G1_r;
                G2_c = G2_r;
            end
            obj@BD_data(G1_r,G2_r,G1_c,G2_c,n);
            obj.XYruns_r = obj.druns('r');
            obj.XYruns_c = obj.druns('c');
            obj.G1_pattern_r = obj.iniPattern(1,'r');
            obj.G2_pattern_r = obj.iniPattern(2,'r');
            obj.G1_pattern_c = obj.iniPattern(1,'c');
            obj.G2_pattern_c = obj.iniPattern(2,'c');
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
            %the lift to the group of either lower- or upper-triangular
            %matrices
           type = strcmp(class(x),'sym');
           if type
              d = poisson_bracket.sym_eye(obj.n).*x;
           else
              d = eye(obj.n).*x;
           end 
           y = obj.Tr(x-d)+d;
        end
        
        function y = Tr(obj,x)
           y = obj.T(x,'r'); 
        end
        function y = Tc(obj,x)
           y = obj.T(x,'c'); 
        end
        function y = Tinvr(obj,x)
           y = obj.Tinv(x,'r'); 
        end    
        function y = Tinvc(obj,x)
           y = obj.Tinv(x,'c'); 
        end
        
        function y = T(obj,x,rc)
            %if rc == 'r', this is for gamma^r
            %if rc == 'c', works for gamma^c
            y = obj.fullT(x,rc); 
        end

        function y = Tinv(obj,x,rc)
            %if rc == 'r', this is for gamma^r
            %if rc == 'c', works for gamma^c
            y = obj.fullT(x,rc,'inv');
        end
        
        function y = fullT(obj,x,rc,varargin)
            inv = false; %indicates whether we use T^{-1} or T
           for i =1:length(varargin) 
               if isa(varargin{i},'char')
                   if strcmpi(varargin{i},'inv')
                      inv = true; 
                   end
               end
           end
           %get the "domain" of the map
           XYruns = obj.getRuns(rc);
           y = poisson_bracket.sym_zeros(obj.n);
           if isempty(XYruns)
              return;
           end
           if ~inv
               dom = obj.getPattern(1,rc);
               start = 1; %the starting position for XYruns for images
               start_im = 3;
           elseif inv
               dom = obj.getPattern(2,rc);
               start = 3;
               start_im = 1;
           end
           x = dom.*x; %first we project onto the union of blocks
           len = length(XYruns(:,1));
           for m = 1:len
               run_length = XYruns(m,start+1)-XYruns(m,start)+1;
               for i = 1:run_length
                  for j = 1:run_length
                      %calculate now positions in matrices
                      ix = XYruns(m,start)+i-1;
                      jx = XYruns(m,start)+j-1;
                      
                      iy = XYruns(m,start_im)+i-1;
                      jy = XYruns(m,start_im)+j-1;
                      
                      y(iy,jy) = x(ix,jx);
                  end
               end
           end
        end

        function y = TTr(obj,x)
           y = obj.TT(x,'r'); 
        end
        function y = TTc(obj,x)
           y = obj.TT(x,'c'); 
        end
        function y = TTinvr(obj,x)
           y = obj.TTinv(x,'r'); 
        end
        function y = TTinvc(obj,x)
           y = obj.TTinv(x,'c'); 
        end
        
        function y = fullTT(obj,x,rc,varargin)
            inv = '';
           for i =1:length(varargin) 
               if isa(varargin{i},'char')
                   if strcmpi(varargin{i},'inv')
                      inv = varargin{i};
                   end
               end
           end
            %realizes the function (1-T)^{-1}
            y = x;
            type = isa(x,'sym');
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
            %t = 1;
            A = obj.fullT(x,rc,inv);
            while true
               y = y + A;
               A = obj.fullT(A,rc,inv);
               %t = t + 1;
               if (A == zero)
                   break;
               end
%                if (t > 10000)
%                   fprintf('The number of iterations exceeded.\n');
%                   fprintf('Review the code for computing (1-T)^{-1}.\n');
%                   break;
%                end
            end
        end
        
        function y = TT(obj,x,rc)
            y = obj.fullTT(x,rc);
        end

        function y = TTinv(obj,x,rc)
            %realizes the function (1-T)^{-1}
            y = obj.fullTT(x,rc,'inv');
        end
        
        function result = proj(obj,A,type,varargin)
            %projects onto the Lie algebra determined by BD-data
            %(or upon its complement)
            %needed only for various tests
            
           %type = 1r,2r,1c,2c 
           %varargin{1} = 'c' for projection onto the complement
           
           complement = sym(zeros(obj.n)); %by default, no complement
           for i = 1:length(varargin)
              if isa(varargin{i},'char')
                  if strcmpi(varargin{i},'c')
                      complement = sym(ones(obj.n));
                  end
              end
           end
           projector = obj.getPattern(str2num(type(1)),type(2));
           if complement(1,1) ~= 0
              projector = complement - projector; 
           end
           result = projector.*A;
        end
    end
    methods (Access=private)
        function XYruns = getRuns(obj,rc)
           if strcmp(rc,'r')
               XYruns = obj.XYruns_r;
           elseif strcmpi(rc,'c')
               XYruns = obj.XYruns_c;
           end
        end
        function G_pattern = getPattern(obj,num,rc)
           if strcmpi(rc,'r')
               if num == 1
                   G_pattern = obj.G1_pattern_r;
               elseif num == 2
                   G_pattern = obj.G2_pattern_r;
               end
           elseif strcmpi(rc,'c')
               if num == 1
                   G_pattern = obj.G1_pattern_c;
               elseif num == 2
                   G_pattern = obj.G2_pattern_c;
               end 
           end
        end
        function G_pattern = iniPattern(obj,num, rc)
            %num -- the number of Gi;
            %rc -- row or column
            G_pattern = poisson_bracket.sym_zeros(obj.n);
            XYruns = obj.getRuns(rc);
            if isempty(XYruns) %case of empty BD data
                return;
            end
            if num == 1
               start = 1;
            elseif num == 2
                start = 3;
            end
            len = length(XYruns(:,1));
            for m = 1:len
               for i = XYruns(m,start):XYruns(m,start+1)
                  for j = XYruns(m,start):XYruns(m,start+1)
                      G_pattern(i,j) = 1;
                  end
               end
            end
        end
    end
end

