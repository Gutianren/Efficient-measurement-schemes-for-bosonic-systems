function [l] = DtoB(d,j,Nq)
    % change d-scale number j to binary number l
    lb = ceil(log2(d)); % number of qubits for each qudit
    Vj = ToDinary(d,j,Nq);
    l = ToNum(2^lb, Vj);
end