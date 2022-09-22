function [C] = update_meas(A,B)
    n = length(A);
    C = zeros(1,n);
    for i = 1:n
        if A(i) == 0
            C(i) = sign(B(i));
        else
            C(i) = sign(A(i));
        end
    end
end