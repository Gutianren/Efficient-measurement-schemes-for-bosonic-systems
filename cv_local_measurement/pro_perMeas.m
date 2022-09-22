function [pro] = pro_perMeas(phi, psi, n, d, dx)

    L = length(psi);
    s1 = reshape(psi,L/d,d)';
    p1 = phi(:, (n-1)*d+1:n*d);
    temp = abs(p1 * s1).^2;
    pro = sum(temp,2) * dx;
    %fprintf('sum(pro)=%f\n',sum(pro));
end