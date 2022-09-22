function [e] = GetObservableExpect()
    global observ phi d
    [m,Nq] = size(observ);
    Nq = Nq-1;
    A = GetGGB(d);
    
    e = 0;
    for i = 1:m
        O = 1;
        for j = 2:Nq+1
            if observ(i,j) == 0
                l = d*d;
            else
                l = observ(i,j);
            end
            O = kron(O,A(:,:,l));
        end
        e = e + observ(i,1) * phi' * O * phi;
    end
    
end