function [x_new, V_new] = GaussianUnitaryTransformation(x,V,d,S)
    x_new = S * x + d;
    V_new = S * V * conj(S');
end