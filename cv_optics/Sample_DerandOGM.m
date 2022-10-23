function [v] = Sample_DerandOGM(measure, T)
    % return an estimator of the ground energy
    
    global observ Nq rho_x rho_V
    
    if T <1
        fprintf('Please input a iteration which is larger than 0.\n');
        return;
    end
    
    
    Sample_mu = zeros(size(measure));
    for t = 1 : T
       Sample_mu(t,:) = GaussianMeasurement(rho_x, rho_V, measure(t,:));
    end
    
    v = 0;
    for j = 1:size(observ,1)
        mu_sum = 0;
        count = 0;  
        for t = 1 : T
            P = measure(t,:);
            if IfCommuteP(observ(j,2:Nq+1), P)
                count = count + 1; % record the number of commuted observables.
                mu = 1;
                for k = 1 : Nq
                    if observ(j,k+1) > 0
                        mu = mu * Sample_mu(t,k)^observ(j,k+1);
                    else
                        mu = mu * Sample_mu(t,k)^(-observ(j,k+1));
                    end
                end
                mu_sum = mu_sum + mu;
            end
        end
        if count ~= 0
           v = v+ mu_sum * 1/count * observ(j,1);
        end
    end
    
end
