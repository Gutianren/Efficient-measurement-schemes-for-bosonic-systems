function [observ] = sort_observ(observ0)
    m = size(observ0,1);
    observ = observ0;
    for i = 1:m
        for j = i+1:m
            if abs(observ(i,1)) < abs(observ(j,1))
                t = observ(i,:); observ(i,:) = observ(j,:); observ(j,:) = t;
            end
        end
    end
end