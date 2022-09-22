function [Sample_mu] = SampleResult_Ground(P, Nq, Coord_phi, d)
    % Sample measured state mu of given GGB string P
    % The result must be Nq elements of -1,1(for sym/asym) or others (for
    % diagonal terms)
    % return a Vector with dimension Nq
    for k = 1:Nq

        [pro,e,u] = pro_perMeas(P(k), Coord_phi, Nq-k + 1, d); %%+1
        j = SampleMeas_global(pro);
        Sample_mu(k) = e(j);
        if k == Nq
            break; % we do not need to update the last qubit.
        end
        Coord_phi = update_gstate(u(:,j), Coord_phi, Nq - k + 1, d);
        %temp = floor((-1 *Sample_mu(k) + 3)/2); %f(1) = 1;f(-1) = 2.
        %Coord_phi(:,2) = Coord_phi(:,2)/sqrt(pro_meas(temp)); %phi_new = M * phi/sqrt(pr(M));
        Coord_phi(:,2) = Coord_phi(:,2) / sqrt(pro(j));
    end
    
end