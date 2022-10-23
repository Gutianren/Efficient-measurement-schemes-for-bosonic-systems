function [O] = rndpx(n,k)
% Generate a N-length p-x string with order = k

O = zeros(1,n);
temp = randperm(n-1+k);
j = 1;
count = 0;
for i = 1:n-1+k
    if temp(i) <= k
        count = count + 1;
    else
        O(j) = count;
        j = j+1;
        count = 0;
    end
end
O(n) = count;
sign = randi(2,1,n)*2-3;
O = O .* sign;
end

