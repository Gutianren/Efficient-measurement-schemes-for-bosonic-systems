function [f] = IfCommuteP(A,P)
    % P can be used to measure A
    n = length(A);
    f = true;
    for i = 1:n
        if sign(A(i)) ~= P(i)
            f = false;
            break;
        end
    end
end
