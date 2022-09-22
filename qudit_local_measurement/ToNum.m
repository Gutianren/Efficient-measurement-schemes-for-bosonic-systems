function [num] = ToNum(d, vector)
    % Change (num)_d to (num)_10

    num = 0;
    for k = 1 : length(vector)
        num = num + vector(k) * d^(length(vector) - k);
    end

end
