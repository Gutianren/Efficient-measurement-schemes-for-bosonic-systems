function [ex] = GetExpectation(rho_x, rho_V, observ)
    M = size(observ,1);
    Nq = size(observ,2)-1;
    ex = 0;
    for i = 1:M
        X = [];
        for j = 1:Nq
            if observ(i,j+1)>0
                for k = 1:observ(i,j+1)
                    X = [X,2*j-1];
                end
            else
                for k = 1:-observ(i,j+1)
                    X = [X,2*j];
                end
            end
        end
        ex = ex + observ(i,1) * GaussianMoment(X,rho_x,rho_V);
    end
end
