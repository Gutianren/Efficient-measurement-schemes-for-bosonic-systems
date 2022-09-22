function [inner] = TraceLocalGen(d,j,l,Q1,Q2)
    % Q1,Q2 are GGB strings (0~d^2-1)
    % 1 <= j,l <= d^2 denote the basis state ( should minus 1)
    % of course we also need the dimension d
    
    Nq = length(Q1);
    j = j - 1;
    l = l - 1;
    
    Vector_j = ToDinary(d,j,Nq);
    Vector_l = ToDinary(d,l,Nq);
    psi = eye(d);
    
    A = GetGGB(d);
    
    inner = 1;
    
    for k = 1:Nq
        if Q1(k) == 0
            Q1(k) = d*d;
        end
        if Q2(k) == 0
            Q2(k) = d*d;
        end
        temp1 = abs(Q1(k));
        temp2 = abs(Q2(k));
        inner = inner * psi(:, Vector_j(k)+1)' * Q1(k)/temp1 * A(:,:,temp1) * Q2(k)/temp2 * A(:,:,temp2) * psi(:, Vector_l(k)+1); %%<j_k |mat_k|l_k>
        if inner == 0
            break;
        end
    end
    
end
