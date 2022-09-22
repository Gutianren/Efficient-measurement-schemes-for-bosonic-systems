function [f] = IfCommute(A,B)
    n = length(A);
    f = true;
    for i = 1:n
        if sign(A(i)) * sign(B(i)) == -1
            f = false;
            break;
        end
    end
    
end
