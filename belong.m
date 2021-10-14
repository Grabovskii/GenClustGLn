function t = belong(t,interval)
%checks if t is between the endpoints
    if isempty(interval)
        t = false;
    else
        t = (t >= interval(1))&&(t <= interval(2));
    end
end

