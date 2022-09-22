function [phi,phi_p] = GetGroundState(x,p)
    global m Nq w gs observ d
    
    %A = exp(1i * p' * x) / sqrt(2*pi);
    %dx = x(2)-x(1);
    for i = 1:Nq
        for j = 1:d
            phi(:, (i-1)*d+j) = Oscillator(w(i),j-1,x)';
            phi_p(:,(i-1)*d+j) = Oscillator(1/w(i),j-1,p)' * (1i)^(j-1);
        end
    end
    %phi_p = A * phi * dx;
end