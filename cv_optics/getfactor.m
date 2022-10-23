factorN = [];
for i = 1:16
    V = squeezedN((ones(i)-eye(i))/i);
    factorN = [factorN;max(max(V))];
end
factor = [];
for i = 1:225
    factor = [factor;factorN(x(i))^(kk(i)/2)];
end

