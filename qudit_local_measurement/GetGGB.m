function [A] = GetGGB(d)
    for i = 1:d^2
        A(:,:,i) = zeros(d);
    end
    
    k = 0;
    for i = 1:d
        for j = i+1:d
            k = k+1;
            A(i,j,k) = 1;
            A(j,i,k) = 1;
            A(i,j,k+d*(d-1)/2) = -1i;
            A(j,i,k+d*(d-1)/2) = 1i;
        end
    end
    
    k = d*(d-1);
    for i = 1:d-1
        k = k+1;
        for j = 1:i
            A(j,j,k) = 1;
        end
        A(i+1,i+1,k) = -i;
        A(:,:,k) = A(:,:,k) * sqrt(2/(i*(i+1)));
    end
    
    A(:,:,d*d) = eye(d);
    
end
