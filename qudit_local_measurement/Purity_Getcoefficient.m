function Purity_Getcoefficient()
    global phi d Nq observ coeff
    A = GetGGB(d);
    ex = 0;
    m = size(observ,1);
    coeff = zeros(m,1);
    for j = 1:m
        Qj = 1;
        n = Nq;
        for k = 2:Nq+1
            p = observ(j,k);
            if p == 0
                p = d^2;
                n = n-1;
            end
            Qj = kron(Qj, A(:,:,p));
        end
        W = 2^n * d^(Nq-n);
        coeff(j) = phi' * Qj * phi / W;
        ex = ex + coeff(j)^2 * W;
    end
    if abs(ex-1) > 1e-4
        fprintf('Warning: rho0 is not a pure state!\n');
        fprintf('Purity(rho0)=%f\n',ex);
    end    
    
end
