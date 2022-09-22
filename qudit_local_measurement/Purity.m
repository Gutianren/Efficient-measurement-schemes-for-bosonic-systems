function [purity] = Purity(rho)
    purity = trace(rho*rho);
end
