function [rho_x, rho_V] = TMSV(r)
    
    v = cosh(2*r);
    u = sqrt(v^2-1);
    I = eye(2);
    Z = [1,0;0,-1];
    rho_x = zeros(4,1);
    rho_V = [v*I,u*Z;u*Z,v*I];
end
