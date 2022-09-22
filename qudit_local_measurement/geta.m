function [ad,a] = geta(d)
    ad = zeros(d,d);
    a = zeros(d,d);
    
    for s = 1:d-1
        ad(s,s+1) = sqrt(s);
        a(s+1,s) = sqrt(s);
    end
end