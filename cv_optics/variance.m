function [var] = variance(pr)

    global observ meas
    
    len = length(pr);
    T = 1000;
    [m, n] = size(observ);
    n = n - 1;
    
    %% diagonal variance %% if there exists an observable Q, which has no rela
    var = 0;
    for j = 1 : m % observ(j)
        temp = 0;
        for k = 1 : len % meas(k)
            if IfCommute(observ(j,2:n+1),meas(:,k))
                temp = temp + pr(k);
            end
        end
        if abs(temp) > 1e-25
            var = var + observ(j,1)^2/temp;
        else
            var = var + observ(j,1)^2*T;
        end
    end
end
