function [v] = Purity_Sample_DerandOGM(observ, Nq, m_new, measure, T, rho)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   derandomize OGM alg by choosing #P = T * pr.
%   At the start, we have observables observ, meaurements meas, the number
%   of qubits Nq, the number of observables m, the number of observables
%   len, and probability for each measurements
%   pr. Here we need to transpose measure to let it in coincidence with the
%   measure of Sample_Derandom.m
%   Output: estimate v after performing T times of sampling.

    %global Nq observ
    global phi d

    v=0;

    if T <1
        display('Please input a iteration which is larger than 0.');
        return;
    end

    [psi, lam] = eig(rho); % A = sum_i lam_i|psi_i><psi_i|
    lam = diag(lam);
    
    for t = 1 : T
    %   Sample_mu(t,:) = SampleResult(measure(t,:), Nq, phi); 
       ps = SampleMeas_global(lam);
       Coord_psi = Coord(psi(:,ps));
       Sample_mu(t,:) = SampleResult_Ground(measure(t,:), Nq, Coord_psi, d); 
    end

    A = GetGGB(d);
    % T_minus = 0;
    % acc = 0;
    v = 0;
    for j = 1 : m_new
        mu_sum = 0;
        count = 0;  
        for t = 1 : T
            P = measure(t,:);
            if IfCommute(observ(j,2:Nq+1), P)
                count = count + 1; % record the number of commuted observables.
                mu = 1;
                for k = 1 : Nq
                    if observ(j,k+1) ~= 0
                        mu = mu * Sample_mu(t,k);
                    end
                end
                mu_sum = mu_sum + mu;
            end
        end
        
%         Qj = [1];
%         for k = 2:Nq+1
%             p = observ(j,k);
%             if p == 0
%                 p = d^2;
%             end
%             Qj = kron(Qj, A(:,:,p));
%         end
%          
%         mu_test = phi' * Qj * phi;
%         mu1 = mu_sum / count;
%         if abs(mu_test-mu1) > 0.1
%             observ(j,:);
%             mu_test;
%             mu1;
%         end

        if count ~= 0
           v = v+ mu_sum * 1/count * observ(j,1);
        end
    end

end

