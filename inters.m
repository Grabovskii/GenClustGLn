function ii = inters(ii1,ii2)
%intersects two intervals that are defined via their endpoints
    if isempty(ii1) || isempty(ii2)
        ii = cell(0);
        return;
    else
        ii = [max(ii1(1),ii2(1)),min(ii1(2),ii2(2))];
        if ii(1) > ii(2)
            ii = cell(0);
        end
    end
end

