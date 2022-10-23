function [V] = squeezedN(beta)
    % squeezed operator = \exp{(\sum_ij \beta_ij*ai*aj - \beta_ij*ai'*aj')/2}
    % using the notation in Gaussian quantum information that p = i(a' - a)
    % suppose beta is real

    N = size(beta,1);
    %beta1 = sqrtm(beta * conj(beta));
    beta1 = beta;
    [U,D] = eig(beta1);
    u = U * diag(cosh(diag(D))) * inv(U);
    v = U * diag(sinh(diag(D))) * inv(U); %* inv(beta1) * beta;
    S = zeros(2*N,2*N);
    for i = 1:N
        for j = 1:N
            S(2*i-1,2*j-1) = (u(i,j) + v(i,j)' + v(i,j) + u(i,j)') / 2;
            S(2*i-1,2*j) = 1i * (u(i,j) + v(i,j)' - v(i,j) - u(i,j)') / 2;
            S(2*i,2*j-1) = -1i * (u(i,j) - v(i,j)' + v(i,j) - u(i,j)') / 2;
            S(2*i,2*j) = (u(i,j) - v(i,j)' - v(i,j) + u(i,j)') / 2;
        end
    end
    V = S * transpose(S);
end
