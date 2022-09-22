function [pro, e, u] = pro_perMeas(measure, Coord_phi, Nq, d)
    % Return the probability of each measurement outcome with a one-site
    % single measure
    % Coord_phi(i,:) = [c1,c2], the coefficient of basis c1 is c2. 
    % measure in {0,1,..,d^2-1} denoting the GGB matrixes
    % if measure is sym/asym, return the probability of measuring plus(+1)
    % if measure is diagonal, return the probability of measuring all kinds
    % of eigenvalues
    % Here Nq is the measured position(also the length of Coord_phi due to
    % the truncation in update_gstate)
    
    if measure == 0
        pro = zeros(1,d);
        e = diag(eye(d));
        pro(1) = 1;
        u = eye(d);
        return;
    end
%     if measure <= d*(d-1)/2 % sym matrixes
%         [i,j] = GetIndex(d, measure);
%         u_plus = zeros(d,1);
%         u_plus(i,1) = 1/sqrt(2);
%         u_plus(j,1) = 1/sqrt(2);
%     elseif measure <= d*(d-1)   % asym matrixes
%         [i,j] = GetIndex(d, measure-d*(d-1)/2);
%         u_plus = zeros(d,1);
%         u_plus(i,1) = 1/sqrt(2);
%         u_plus(j,1) = 1i/sqrt(2);
%     else % diagonal matrixes, have all kinds of eigenvectors
%         max_j = measure - d*(d-1) + 1; 
%     end
    A = GetGGB(d);
    M = A(:,:,measure);
    [u, e] = eig(M);
    e = diag(e);
    Coord_len = size(Coord_phi,1); %%number of non-zero baisis
    %%calcuate |sum_{j} alpha_j <u_plus|j>|^2.
    identity = eye(d);
    pro = zeros(d,1);
    
    for l = 1:d
        for j = 1 : Coord_len
            CoordAfirst = floor(Coord_phi(j,1)/d^(Nq-1));
            CoordArest = mod(Coord_phi(j,1), d^(Nq-1)); %%have revised***

            sum_a = Coord_phi(j,2)';
            sum_a = sum_a * identity(:,CoordAfirst+1)' * u(:,l);
            sum_b = 0;
            if sum_a ~= 0
                for k = 1 : Coord_len
                    CoordBfirst = floor(Coord_phi(k,1)/d^(Nq-1));
                    CoordBrest = mod(Coord_phi(k,1), d^(Nq -1));
                    if CoordBrest == CoordArest
                        sum_b = sum_b + Coord_phi(k,2)* u(:,l)' * identity(:,CoordBfirst + 1);
                    end
                end
            end
            
            pro(l) = pro(l) + sum_a * sum_b;
        end
    end
    
end