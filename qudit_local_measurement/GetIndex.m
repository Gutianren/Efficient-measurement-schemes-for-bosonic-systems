function [i,j] = GetIndex(d, m)
    % Get the index of symmetric GGB matrix
    if m > 0 && m <= d*(d-1)
       i = 1;
       while m > d-i
           m = m - (d-i);
           i = i + 1;
       end
       j = m + i;
       return;
    else
       fprintf('Invalid input for GetIndex');
       pause;
    end
end