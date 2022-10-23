function [result, x_new, V_new] = SingleGaussianMeasurement(x, V, meas)
    % Measuring q or p for the first mode of the Gaussian state with heterodyne detection
    % meas = 0 or 1 denoting measuring on q or p
    N = size(x,1) / 2;
    A = V(1:2,1:2);
    B = V(3:end,3:end);
    C = V(1:2,3:end);
    xA = x(1:2);
    xB = x(3:end);
    if meas >= 0
        result = mvnrnd(x(1),V(1,1));
        xm = [result;0];
        P = diag([1,0]);
    else
        result = mvnrnd(x(2),V(2,2));
        xm = [0;result];
        P = diag([0,1]);
    end
    
    PAPinv = pinv(P*A*P);
    x_new = xB - C' * PAPinv * (xA - xm);
    V_new = B - C' * PAPinv * C;
    
end
