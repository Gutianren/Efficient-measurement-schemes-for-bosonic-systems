function [newstate] = update_gstate(j_x, phi, psi, n, d)
    L = length(psi);
    psi = reshape(psi,L/d,d);
    p1 = phi(j_x, (n-1)*d+1:n*d)';
    newstate = psi * p1;
    newstate = newstate / norm(newstate);
end