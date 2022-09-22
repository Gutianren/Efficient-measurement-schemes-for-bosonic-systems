function [coe, idx] = decompose(d,A)
    % decompose a Hermitian matrix into Lambda-Matrixes
    % indexes:
    % 1~d(d-1)/2 = Lambda_sym
    % d(d-1)/2+1 ~ d(d-1) = Lambda_asym, 
    % d(d-1)+1 ~ d^2-1 = Lambda_l
    % d^2 = I
    pos = 0;
    for i = 1:d
        for j = i+1:d
            pos = pos + 1;
            coes(pos) = real(A(i,j));
            coes(pos+d*(d-1)/2) = -imag(A(i,j));
        end
    end
    
    tr(1) = A(1,1);
    
    for i = 2:d
        tr(i) = tr(i-1) + A(i,i);
    end
    
    for l = 1:d-1
        coes(d*(d-1)+l) = (tr(l)-l*A(l+1,l+1)) / sqrt(2*l*(l+1));
    end
    
    coes(d*d) = tr(d)/d;
    
    L = 0;
    for i = 1:d*d
        if coes(i) ~= 0
            L = L + 1;
            coe(L) = coes(i);
            idx(L) = i;
        end
    end
    
    coe = coe';
    idx = idx';
end