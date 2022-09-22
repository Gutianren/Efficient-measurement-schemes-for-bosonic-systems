function [gs_new] = x_shift(gs, a, q)
    % Shifting the state with exp(ip(a*dx)) on the q-th qubit
    % in the energy space
    % a must be a slight shift(<< xmax / dx)
    global d Nq phi dx
    
    coef = zeros(d,d); % coeff(n,m) denote sum_x( phi_n(x) * phi_m(x+a)' )
    phi_a = [phi(a+1:end,:); phi(1:a,:)];
    sum(abs(phi_a).^2)
    for n = 1:d
        for m = 1:d
            coef(n,m) = phi_a(:,(q-1)*d+m)' * phi(:,(q-1)*d+n) * dx;
        end
    end
    coef * coef'
    gs_new = zeros(d^Nq,1);
    for j = 1:d^Nq
        vj = ToDinary(d,j-1,Nq);
        for k = 1:d^Nq
            vk = ToDinary(d,k-1,Nq);
            f = find(vj ~= vk);
            if isempty(f) || length(f) == 1 && f == q
                gs_new(j) = gs_new(j) + gs(k) * coef(vj(q)+1, vk(q)+1);
            end
        end
    end
    sum(abs(gs_new).^2)
end