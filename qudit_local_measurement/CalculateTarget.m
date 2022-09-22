function [target] = CalculateTarget(d,j, Q1, Q2)
    % calculate the state when Q act on it
    Nq = length(Q1);
    j = j - 1;
    Vector_j = ToDinary(d,j,Nq);
    psi = eye(d);
    Vector_l = zeros(1,Nq);
    A = GetGGB(d);
    for i = 1:Nq
        if Q1(i) == 0
            Q1(i) = d^2;
        end
        if Q2(i) == 0
            Q2(i) = d^2;
        end
        
        p1 = A(:,:,Q1(i)) * A(:,:,Q2(i)) * psi(:,Vector_j(i)+1);
        if ~isempty(find(p1~=0, 1))
            Vector_l(i) = find(p1~=0,1)-1;
        else
            target = -1;
            return;
        end
    end
    target = ToNum(d,Vector_l) + 1;
end