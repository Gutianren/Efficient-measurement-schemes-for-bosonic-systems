function [result] = GaussianMeasurement(x,V,meas)
    % Measuring all N quantities of Gaussian state (x,V)
    % meas = 0 or 1 denoting x or p measurement
    
    N = size(x,1) / 2;
    result = zeros(N,1);
    for i = 1:N
        if meas(i) ~= 0
            [result(i), x, V] = SingleGaussianMeasurement(x,V,meas(i));
        else
            [result(i), x, V] = SingleGaussianMeasurement(x,V,1);
        end
    end
    
end