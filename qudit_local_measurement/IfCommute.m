function f = IfCommute(A,B)
    n = length(A);
    f = true;
    for i = 1:n
        if A(i) ~= B(i) && A(i) ~= 0 && B(i) ~= 0
            f = false;
            break;
        end
    end
end