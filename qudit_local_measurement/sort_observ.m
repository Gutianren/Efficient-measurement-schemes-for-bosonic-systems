function [coe,observ] = sort_observ(coe0,observ0)
    m = length(coe0);
    coe = coe0;
    observ = observ0;
    for i = 1:m
        for j = i+1:m
            if abs(coe(i)) < abs(coe(j))
                t = coe(i); coe(i) = coe(j); coe(j) = t;
                t = observ(i,:); observ(i,:) = observ(j,:); observ(j,:) = t;
            end
        end
    end
end