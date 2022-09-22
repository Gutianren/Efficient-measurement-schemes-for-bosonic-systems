Repeat = 10;

TT = [100,200,500,1000,2000,5000];
%T = 10000;

for d = 3:7
    for t = 1:length(TT)
        T = TT(t);
        main(d,T,Repeat);
    end
end
