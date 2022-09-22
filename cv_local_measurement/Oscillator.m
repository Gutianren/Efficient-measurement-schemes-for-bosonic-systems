function [phi] = Oscillator(w,n,x)
    % Get the n-th level of Oscillator with frequency w
    % In the x space
    phi = sqrt(1/(2^n*factorial(n))) * (w/pi)^(1/4) * exp(-w*x.^2/2) .* hermiteH(n, sqrt(w) * x);
end
