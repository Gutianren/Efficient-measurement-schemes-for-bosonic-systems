function Purity_Getobserv()

    global Nq observ d
    m = d^(2*Nq);
    observ = zeros(m, Nq+1);
    %observ(:,1) = ones(m,1) / 2^Nq;
    for j = 1:m
        j1 = j-1;
        for k = 1:Nq
            observ(j,Nq+2-k) = mod(j1,d^2);
            j1 = floor(j1/d^2);
        end
        n = sum(abs(observ(j,2:Nq+1))>1e-10);
    end

end
