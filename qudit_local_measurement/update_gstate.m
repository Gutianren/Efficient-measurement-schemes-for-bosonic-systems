function [ Coord_new ] = update_gstate( u_sample, Coord_phi, Nq, d )
    

    identity = eye(d);

    Coord_len = size(Coord_phi,1); %%number of non-zero baisis
    
    % calcuate |sum_{j} alpha_j <u_plus|j>|^2.
    cur = 0;

    for j = 1 : Coord_len
        alp_j = Coord_phi(j,2);
        firstQ = floor(Coord_phi(j,1)/d^(Nq-1));
        amplitude = alp_j * u_sample' * identity(:,firstQ+1);
     %  fprintf('%d: alp_%d = %f%+fj, firstQ is %d\n', j, Coord_phi(j,1), real(alp_j),imag(alp_j), firstQ);
        if amplitude ~= 0
            cur = cur + 1;
            basis_new = mod(Coord_phi(j,1),d^(Nq-1));
            Coord_new(cur,1) = basis_new;
            Coord_new(cur,2) = amplitude;
        end
    end
    
    % Clear the same terms
    len_new = cur;
    cur = 0;
    for j = 1 : len_new
        if Coord_new(j,1) == -1
            continue;
        end
        cur = cur + 1;
        Coord_phi(cur,1) = Coord_new(j,1);
        Coord_phi(cur,2) = Coord_new(j,2);
        for k = j+1 : len_new   
            if Coord_new(k, 1) == Coord_new(j,1)
                Coord_new(k, 1) = -1; %lable this basis has been visited.
                Coord_phi(cur, 2) = Coord_phi(cur,2) + Coord_new(k,2);
            end
        end
    end
    Coord_new = Coord_phi(1:cur,:);

end