function [result] = GaussianMeasurement2(x, V, meas)
    %% It is a test program
    % Measuring N variables(q or p) of the Gaussian state with heterodyne detection
    % meas = 0 or 1 denoting measuring on q or p
    
    N = size(x,1) / 2;
    idx = (0:N-1)*2 + meas + 1;
    x_sub = x(idx);
    V_sub = V(idx,idx);
    result = mvnrnd(x_sub, V_sub);

end
