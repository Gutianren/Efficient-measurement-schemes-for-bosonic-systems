function [ex] = GaussianMoment(X, mu, Cov)
%   Calculate the moment of gaussians
%   X are the labels of the multiplied variables
%   mu is the expectation
%   Cov is the covariance matrix
%   Using X_i = Y_i + mu_i to calculate, where Y_i are the central
%   variables

    L = length(X);
    ex = 0;
    for T = 1:2^L
        t = dec2bin(T-1,L);
        tempex = 1;
        temparr = [];
        for i = 1:L
            if t(i) == '0'  % choose Y_i = X_i - mu_i
                temparr = [temparr, X(i)];
            else    % choose mu_i
                tempex = tempex * mu(X(i));
            end
        end
        ex = ex + tempex * GaussianCentralMoment(temparr, Cov);
    end

end