function [psi] = BinaryGroundState(d,phi)
    % encode the ground state into qubits (compact encoding)
    % phi is a d^Nq array denoting the coefficients
    
    Nq = floor(log(length(phi)) / log(d));
    lb = ceil(log2(d));
    psi = zeros(2^(Nq*lb),1);
    for j = 1:length(phi)
        if phi(j) ~= 0
            l = DtoB(d,j-1,Nq) + 1;
            psi(l) = phi(j);
        end
    end
end