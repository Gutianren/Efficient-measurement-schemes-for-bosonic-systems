function measP = sort_meas(mp)
    % sort the measurements according to the absolute value of the
    % probability
    [n,m] = size(mp);
    n = n - 1;
    measP = mp;
    
    for i = 1:m
        for j = i+1:m
            if abs(measP(n+1,i)) < abs(measP(n+1,j))
                t = measP(:,i); measP(:,i) = measP(:,j); measP(:,j) = t;
            end
        end
    end
    
end