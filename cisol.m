classdef cisol
    %this class realizes the hidden functions denoted as c
    %WARNING: the indexation of the cell array C is shifted
    %with respect to the Misha's indexation. C{1} corresponds to detX
    %and C{n+1} to detY. There are functions get_c and c with the correct
    %indexation; c(0) = detX and c(n) = detY
    
    properties
        C
        n
    end
    
    methods
        function obj = cisol(n)
            %function produces isolated functions c
            %NOTE: the indexation is shifted 
            X = sym('X',[n n]);
            Y = sym('Y',[n n]);
            t = sym('t');
            P = det(X+t*Y);
            obj.n = n;
            obj.C = cell(1,n+1);
            obj.C{1} = subs(P,t,sym('0'));
            for k = 2:(n+1)
                P = diff(P,t);
                obj.C{k} = subs(((-1)^((k-1)*(n-1)))*(1/factorial(k-1))*P,t,sym('0'));
            end
        end
        function c = get_c(obj,k)
            %this makes indexation as in article; i.e., starting from 0
            if (k > -1)&&(k < obj.n+2)
                c = obj.C{k+1};
            else
                c = sym('1');
            end
        end
        function cc = c(obj,k)
            cc = obj.get_c(k);
        end
    end
end

