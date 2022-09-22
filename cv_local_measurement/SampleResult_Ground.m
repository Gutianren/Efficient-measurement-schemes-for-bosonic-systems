function [Sample_mu] = SampleResult_Ground(x, p, P, phi, phi_p, psi, d, Nq)
    % Sample one result from measurement P
    % x,p denoting the sample basis
    % phi and phi_p are connected through Fourier Transformation
    % psi is the ground state in the energy basis
    
    for k = 1:Nq
        if P(k) == 1
            pro = pro_perMeas(phi, psi, k, d, x(2)-x(1));
            j = SampleMeas_global(pro);
            Sample_mu(k) = x(j);
            if k == Nq
                break;
            end
            psi = update_gstate(j, phi, psi, k, d);
        else
            pro = pro_perMeas(phi_p, psi, k, d, p(2)-p(1));
            j = SampleMeas_global(pro);
            Sample_mu(k) = p(j);
            if k == Nq
                break;
            end
            psi = update_gstate(j, phi_p, psi, k, d);
        end
    end
    
end