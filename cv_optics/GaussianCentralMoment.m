function [ex] = GaussianCentralMoment(X, Cov)
%   Calculate the central moment of gaussians
%   X are the labels of the multiplied variables
%   Cov is the covariance matrix
    L = size(X,2);
    ex = 0;
    if mod(L,2) == 1
        return;
    end
    if L == 0
        ex = 1;
        return;
    end
    for i = 2:L
        sub = GaussianCentralMoment([X(2:i-1),X(i+1:L)],Cov);
        ex = ex + Cov(X(1),X(i)) * sub;
    end
end

