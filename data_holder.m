classdef data_holder < handle
    %this class is designed to perform various manipulations 
    %with the data on the app level
    %unfortunately, the performance of the app designer is so low
    %that one has to write as much code aside as possible
    
    %inherited handle in order to access it in Poisson Calculator
    
    properties
        clust_str
        %clust_str_r
        only_quiver
        %U
        %U0
        e %error handler, which contains all instances of possible errors
        birat %holder for (W,Z)-structure and U
    end
    
    methods
        function obj = data_holder(varargin)
            len = length(varargin);
            %default
           build_phi = 0;
           build_c = 0;
           bracket = 1;
           build_ghf = 1;
           build_birat = 0;
           obj.only_quiver = 0;
            
            for i = 1:len
               if ischar(varargin{i})
                   if strcmp(varargin{i},'G1_r')
                       G1_r = error_handler.extract_G(varargin{i+1});
                   elseif strcmp(varargin{i},'G2_r')
                       G2_r = error_handler.extract_G(varargin{i+1});
                   elseif strcmp(varargin{i},'G1_c')
                       G1_c = error_handler.extract_G(varargin{i+1});
                   elseif strcmp(varargin{i},'G2_c')
                       G2_c = error_handler.extract_G(varargin{i+1});
                   elseif strcmp(varargin{i},'build_phi')
                       build_phi = varargin{i+1};
                   elseif strcmp(varargin{i},'build_c')
                       build_c = varargin{i+1};
                   elseif strcmp(varargin{i},'only_quiver')
                       obj.only_quiver = varargin{i+1};
                   elseif strcmp(varargin{i},'bracket')
                       bracket = varargin{i+1};
                   elseif strcmp(varargin{i},'n')
                       n = varargin{i+1};
                   elseif strcmp(varargin{i},'build_birat')
                       build_birat = varargin{i+1};
                   elseif strcmp(varargin{i},'type')
                       type = varargin{i+1};
                   end
               end
            end
            if obj.only_quiver
               build_phi = 0;
               build_c = 0;
               bracket = 0;
               build_ghf = 0;
               build_birat = 0;
            else
                build_ghf = 1;
            end
            if true %check the type
               obj.clust_str = double_full(G1_r,G2_r,G1_c,G2_c,n,...
                   'build_phi',build_phi,'build_c',build_c,...
                   'build_ghf',build_ghf,'bracket',bracket);
               if G1_r(1) ~= 0 && build_birat == 1
                    obj.birat = testbirat(G1_r,G2_r,G1_c,G2_c,n,'double',obj.clust_str);
                    %build_birat = 0; <- как это сюда попало? зачем это?
               end
            end
            obj.e = error_handler('n',n,'build_ghf',build_ghf,'build_c',build_c,...
                'bracket',bracket,'build_p',build_phi,'build_birat',build_birat); %a handler of some errors
        end
        
        function output = show_U(obj)
           if isempty(obj.birat)
               output = {'The birational map has not been constructed.'};
               return;
           else
               n = obj.clust_str.q.n;
               output = cell(n + 2,1);
               output{1} = 'U = ';
               output(2:n+1) = obj.disp2cell(obj.birat.U);
               output{n+2} = ' ';
           end
        end
        
        function output = show_invU(obj)
           if isempty(obj.birat)
               output = {'The birational map has not been constructed.'};
               return;
           else
               if isempty(obj.birat.U_1)
                   obj.birat.invU();
               end
               n = obj.clust_str.q.n;
               output = cell(n + 2,1);
               output{1} = 'inv(U) = ';
               output(2:n+1) = obj.disp2cell(obj.birat.U_1);
               output{n+2} = ' ';
           end
        end
        
        function output = show_U0(obj)
           if isempty(obj.birat)
               output = {'The birational map has not been constructed.'};
               return;
           else
               n = obj.clust_str.q.n;
               output = cell(n + 2,1);
               output{1} = 'U0 = ';
               output(2:n+1) = obj.disp2cell(obj.birat.U0);
               output{n+2} = ' ';
           end
        end
        
        function output = display_L(obj)
           if obj.clust_str.build_ghf_ == 0
               output = cell(1);
               output{1} = 'The functions have not been constructed. Error!';
               return;
           end
           count = length(obj.clust_str.s.L);
           Llen = zeros(count,1);
           output_length = 0;
           for i = 1:count
               [Llen(i),~] = size(obj.clust_str.s.L{i});
               output_length = output_length + Llen(i) + 2;
           end
           output = cell(output_length, 1);
           output{1} = 'L(1) = ';
           output(2:Llen(1)+1) = obj.disp2cell(obj.clust_str.s.L{1});
           output{Llen(1)+2} = ' ';
           new_position = Llen(1) + 3;
           for i = 2:count
              output{new_position} = append('L(',num2str(i),') = ');
              output(new_position+1:Llen(i)+new_position) = obj.disp2cell(obj.clust_str.s.L{i});
              output{new_position+Llen(i)+1} = ' ';
              new_position = new_position+Llen(i)+2;
           end
        end
        
        function f = get_f_WZ(obj,l,i,varargin)
           %all errors should be caught before calling this function
            j = cell(0);
            if ~isempty(varargin)
                j = varargin{1};
            end
           if  strcmpi('c',l)
                f = obj.birat.bd_r.c.get_c(i);
           elseif contains('gfhp',l,'IgnoreCase',true)
               f = obj.birat.bd_r.get_f(l,i,j);
           elseif strcmpi('w',l)
               f = obj.birat.bd_r.b.X(i,j);
           elseif strcmpi('z',l)
               f = obj.birat.bd_r.b.Y(i,j);
           elseif strcmpi('u',l)
               f = obj.birat.U(i,j);
           elseif strcmpi('u0',l)
               f = obj.birat.U0(i,j);
           elseif strcmpi(l,'a')
               %f = obj.birat.U0(i,j)*obj.birat.bd_r.s.g(obj.birat.root_r+1,1);
               %f = simplifyFraction(f,'Expand',true);
               [f,~] = numden(obj.birat.U0(i,j)); 
               %I'm not sure what works faster and am reluctant to check
           end
        end
        
        function f = show_get(obj,l,i,varargin)
           %all errors should be caught before calling this function
           %it's called in show_fn
            j = cell(0);
            if ~isempty(varargin)
                j = varargin{1};
            end
           if contains('gfhp',l,'IgnoreCase',true)
               f = obj.clust_str.get_m(l,i,j);
           elseif strcmpi('u',l)
               f = obj.birat.U(i,j);
           elseif strcmpi('u0',l)
               f = obj.birat.U0(i,j);
           elseif strcmpi('c',l)
               f = obj.clust_str.c.get_c(i);
           elseif strcmpi('a',l)
               f = obj.get_f_WZ('a',i,j); 
           end
        end
        
        function output = show_fn(obj,str)
            %there should be a more clever way of catching u0
            %but I'm reluctant to search for it
            if contains(str,'u0','IgnoreCase',true)
               txt = {'u0'};
               str = erase(str,'u0');
            else 
                txt = extract(str,lettersPattern);
            end
            inds = extract(str, digitsPattern);
            inds = cellfun(@str2num,inds); 
            
            if length(inds) == 1
                output = obj.e.check('show',txt{1},inds(1));
            elseif length(inds) == 2
                output = obj.e.check('show',txt{1},inds(1),inds(2));
            else
                output = {'The type of variable is not indicated.'};
            end
            
            if ~isempty(output)
                return;
            end
            
            if ~strcmp(txt{1},'c')
                f = obj.show_get(txt{1},inds(1),inds(2));
                [n,~] = size(f);
                output = cell(n+2,1);
                output{1} = append(txt{1},'(',num2str(inds(1)),',',num2str(inds(2)),') = ');
                output(2:n+1) = obj.disp2cell(f);
                output{n+2} = ' ';
            else
                output = cell(3,1);
                output{1} = append(txt{1},'(',num2str(inds(1)),') = ');
                output(2) = obj.disp2cell(obj.show_get('c',inds(1)));
                output{3} = ' ';
            end
        end
        
        function f = form(obj,f)
           [N,~] = size(f);
           for i = 1:N
               for j = 1:N
                  if f(i,j) ~= 0
                      f(i,j) = 1;
                  end
               end
           end
        end
        
        function f = XY2WZ(obj,f)
            if isempty(obj.birat)
                f = cell(0);
                return;
            end
            f = subs(f,[obj.clust_str.s.X,obj.clust_str.s.Y],[obj.birat.UW,obj.birat.UZ]);
        end
        
        function f = WZ2XY(obj,f)
            if isempty(obj.birat)
                f = cell(0);
                return;
            end 
            f = subs(f,[obj.birat.bd_r.b.X,obj.birat.bd_r.b.Y],...
                [obj.birat.U_1X,obj.birat.U_1Y]);
        end
        
        function output = bracket(obj,str1,str2,varargin)
           %this method is mainly for data_holder object if
           % the app version is run
           [lets,inds] = obj.read_off_bracket(str1,str2);
           
           for i =1:length(lets)
               if isempty(lets(i))
                   output{1} = 'The type of one of the variables has not been specified.';
                   return;
               end
           end
           
           if isempty(obj.clust_str.b)
              output{1} = 'The bracket has not been constructed.';
              return;
           end
           %let's determine whether the bracket is taken in (X,Y) or (W,Z)
           bracket_type = 'bracket';
           len = length(varargin);
           for i = 1:len
               if strcmpi(varargin{i},'wz')
                   bracket_type = 'wz';
                   break;
               end
           end
           
           output = obj.e.check(bracket_type,lets{1},inds{1},inds{2});
           if ~isempty(output)
               return;
           else
               output = obj.e.check(bracket_type,lets{2},inds{3},inds{4});
               if ~isempty(output)
                   return;
               end
           end
           
           log_situation = 0;
           U_bracket = 0;
           bracket_U_U = 0;
           if ~strcmpi(bracket_type,'wz')
               f1 = obj.get_f(lets{1},inds{1},inds{2});
               f2 = obj.get_f(lets{2},inds{3},inds{4});
           else
               f1 = obj.get_f_WZ(lets{1},inds{1},inds{2});
               f2 = obj.get_f_WZ(lets{2},inds{3},inds{4});
           end
           for i = 1:len
              if strcmpi(varargin{i},'log')
                 f1 = log(f1);
                 f2 = log(f2);
                 log_situation = 1;
              elseif strcmpi(varargin{i},'U_bracket')
                 U_bracket = 1; %apply U after the result is computed
                 if isempty(obj.birat)
                    output{1} = 'The birational map has not been constructed.';
                    return;
                 end
              elseif strcmpi(varargin{i},'bracket_U_U')
                 bracket_U_U = 1; %apply U after the result is computed
                 if isempty(obj.birat)
                    output{1} = 'The birational map has not been constructed.';
                    return;
                 end 
              end
           end
          %separation of indices of lets ~= 'c'
           sp = cell(2,1);
           for k = 1:2
               if strcmp(lets{k},'c')
                   sp{k} = '';
               else
                   sp{k} = ',';
               end
           end
           index1 = append('(',num2str(inds{1}),sp{1},num2str(inds{2}),')');
           index2 = append('(',num2str(inds{3}),sp{2},num2str(inds{4}),')');
           
           um = cell(2,1);
           if U_bracket
               um{1} = 'U';
               um{2} = '';
           else
              um{1} = '';
              um{2} = '';
           end
           
           muu = cell(2,1);
           if bracket_U_U
               muu{1} = 'U(';
               muu{2} = ')';
           else
              muu{1} = '';
              muu{2} = '';
           end               
           
           lm = cell(2,1);
           if log_situation
              lm{1} = 'log(';
              lm{2} = ')';
           else
               lm{1} = '';
               lm{2} = '';
           end
           
           

           msgbr = append(um{1},'{',muu{1},lm{1},lets{1},index1,lm{2},muu{2},', ', muu{1},lm{1},lets{2},index2,lm{2},muu{2},'} = '); 
           
           if strcmpi(bracket_type,'wz')
               result = obj.birat.bd_r.b.bracket(f1,f2);
           elseif ~bracket_U_U
               result = obj.clust_str.b.bracket(f1,f2);
           else
              result = obj.birat.bd_r.b.bracket(obj.birat.apply_U(f1),obj.birat.apply_U(f2));
           end
           if U_bracket
              result = obj.birat.apply_U(result);
           end
           result = simplifyFraction(result,'Expand',true);
           msgbr = append(msgbr, char(result));
           output = cell(2,1);
           output{1} = msgbr;
           output{2} = ' ';
        end

        function output = show_R(obj,type,varargin)
            case_wz = 0;
            case_diff = 0;
            form = false;
            for i =1:length(varargin)
                if ischar(varargin{i})
                    if strcmpi(varargin{i},'wz')
                        case_wz = 1;
                    end
                    if strncmpi(varargin{i},'diff',4)
                        case_diff = 1;
                    end
                end
                if isa(varargin{i},'logical')
                   form = varargin{i}; 
                end
            end
           if isempty(obj.clust_str.b)
               output = cell(1,1);
               output{1} = 'The bracket has not been constructed.';
               return;
           elseif (case_wz||case_diff) && ~obj.e.build_birat
               output = cell(1,1);
               output{1} = 'The (W,Z)-variety has not been constructed.';
               return;
           end
           A = sym('A',obj.clust_str.q.n);
           if strncmpi(type,'row',3)
               if ~case_wz && ~case_diff
                   res = obj.clust_str.b.Rr(A);
                   msg = append('R_r(A) = ');
               elseif case_diff
                   res = obj.clust_str.b.Rr(A) - obj.birat.bd_r.b.Rr(A);
                   msg = append('R_r,xy(A) - R_r,wz(A) = ');
               else
                   res = obj.birat.bd_r.b.Rr(A);
                   msg = append('R_r(A) = ');
               end
               
           elseif strncmpi(type,'col',3)
               if ~case_wz
                   res = obj.clust_str.b.Rc(A);
               else
                   res = obj.birat.bd_r.b.Rc(A);
               end
               msg = append('R_c(A) = ');
           else
               output = cell(1); 
               output{1} = 'Check out the code. ERROR.';
               return;
           end
           if form
              res = obj.form(res); 
           end
           n = obj.clust_str.q.n;
           output = cell(n+2,1);
           output{1} = msg;
           output(2:n+1) = obj.disp2cell(res);
           output{n+2} = ' ';
        end
        
        function output = mutate_seq(obj,str)
            %mutates along the sequence prescribed in the string str
            %forms a string _output_ that will be then sent to the 
            %output text area
            %errormsg = 'Error in the syntax. No mutation produced.';
            str = erase(str,'hi'); %if 'phi' is entered, change into 'p'
            lets = extract(str,lettersPattern);
            inds = extract(str,digitsPattern);
            len = length(lets);
            %if 2*len ~= length(inds)
            %    output{1} = errormsg;
            %    return;
            %end
            output = cell(2*len,1);
            inds = cellfun(@str2num,inds);
            %you can see many commented 2*i-1. this is because I changed
            %the output in the relevant textbox upside down
            for i = 1:len
                if obj.only_quiver
                    check = obj.e.check('quiver',lets{i},inds(2*i-1),inds(2*i));
                else
                   check = obj.e.check('mutation',lets{i},inds(2*i-1),inds(2*i));
                end
                if ~isempty(check)
                    %output(2*i-1) = check;
                    output(2*i) = check;
                    continue;
                end
               index = append('(',num2str(inds(2*i-1)),',',num2str(inds(2*i)),')');
               if obj.clust_str.q.isfrozen(obj.clust_str.q.ind(inds(2*i-1),inds(2*i),lets{i}))
                  %output{2*i-1} = append('The variable ',lets{i},index,' is frozen!');
                  output{2*i} = append('The variable ',lets{i},index,' is frozen!');
                  %output{2*i} = ' ';
                  continue;
               end
               %isempty(output{2*i-1}) -- by this I mean that no mistake
               %has been caught (therefore, the input is mutable)
               if obj.only_quiver && isempty(output{2*i-1})
                   obj.clust_str.q.mutate(lets{i},inds(2*i-1),inds(2*i));
               elseif ~obj.only_quiver && isempty(output{2*i-1})
                   obj.clust_str.mutate(lets{i},inds(2*i-1),inds(2*i));
                   %mut = append('Mutated at ',lets{i},index,';\n');
                   f = char(obj.clust_str.get_f(lets{i},inds(2*i-1),inds(2*i)));
                   %output = append(mut,f,'\n');
                   %output{2*i-1} = append(lets{i},'''',index,'= ',f);
                   %output{2*i} = ' ';
                   output{2*i} = append(lets{i},'''',index,'= ',f);
                   output{2*i-1} = ' ';
               end
            end
            output = output(~cellfun('isempty',output));
            %in case we wanted only quiver mutations, correct the output
            if obj.only_quiver
                %output = [output; {'The quiver mutations have been performed. Plot the quiver to see the new seed.'}];%; {' '}];
                output = [{'The quiver mutations have been performed. Plot the quiver to see the new seed.'}; output];
            end
            
            output = flip(output); %just because of my new convention
        end
        function f = get_f(obj,type,i,varargin)
            %get a function to compute the bracket in (X,Y)
            if ~isempty(varargin)
                j = varargin{1};
            else
                j = cell(0);
            end
            if isempty(j)
                f = obj.clust_str.c.get_c(i);
            else
               if strcmpi(type,'g')
                  f = obj.clust_str.s.g(i,j); 
               elseif strcmpi(type,'f')
                  f = obj.clust_str.d.f(i,j);
               elseif strcmpi(type,'h')
                  f = obj.clust_str.s.h(i,j);
               elseif strcmpi(type,'phi')||strcmpi(type,'p')
                   f = obj.clust_str.d.phi(i,j);
               elseif strcmpi(type,'x')
                   f = obj.clust_str.s.X(i,j);
               elseif strcmpi(type,'y')
                   f = obj.clust_str.s.Y(i,j);
               end
            end
        end
        function output = disp2cell(obj,f)
            %an expression or a matrix f is transformed to a cell
            %array of chars suitable for outputtext in appdesigner
            if isempty(f)
                output = cell(1);
                output{1} = char();
                return;
            end
            displayed = evalc('disp(f)');
            displayed = splitlines(displayed); 
            output = displayed(1:length(displayed)-2);
        end
    end
    methods( Access= private)
        function [lets, inds] = read_off_bracket(obj,str1,str2)
            
           %converts a string into an understandable array of indices 
           %for the bracket
           %it could be done neat, but, alas, we have c-functions and u0
           %that spoil the laconic code.
           
           %SO, this function always returns a cell array lets of size 2
           %with possibly empty entries (then we throw an error)
           %AND it returns an array inds of size 4, regardless of whether
           %we use c-functions or not
           
           str = cell(2,1); 
           str{1} = lower(str1); str{2} = lower(str2);
           lets = cell(2,1);
           inds = cell(4,1);
           for i =1:2
               if contains(str{i},'u0')
                   %erase so that in digitsPattern we don't catch extra 0
                   str{i} = erase(str{i},'u0');
                   lets(i) = {'u0'};
               else
                   %just in case the user left garbage letters in the input
                   bb = extract(str{i},lettersPattern);
                   if ~isempty(bb)
                       lets(i) = bb(1);
                   end
               end
               indsx = extract(str{i},digitsPattern);
               if ~isempty(indsx)
                   inds(2*i-1:length(indsx)+ 2*i-2) = indsx;
               end
           end
           
           for i = 1:4
              if ~isempty(inds{i})
                  inds{i} = str2num(inds{i});
              end
           end
           
        end
    end
end

