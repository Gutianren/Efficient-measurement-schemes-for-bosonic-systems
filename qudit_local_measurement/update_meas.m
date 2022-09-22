function C = update_meas(A,B)
    n = length(A);
    for i = 1:n
        if A(i) == 0
            C(i) = B(i);
        else
            C(i) = A(i);
        end
    end
end