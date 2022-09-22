function GetObserv(file)
    
    global Nq observ w m V
    % x : 1, 2, ... , p : -1, -2, -3, ...
    % I : 0
    MaxNq = 20;
    
    Nq = 0;
    file = strcat(file,'.txt');
    f1 = fopen(file);
    j = 0;
    observ = zeros(1,MaxNq);
    while ~feof(f1)
        j = j+1;
        fline = fgetl(f1);
        datal = str2num(fline);
        L = length(datal);
        observ(j,1) = datal(L);
        
        if L-1 == 2 && datal(1) == datal(2)
           i = datal(1);
           w(i) = sqrt(2*datal(L));
        end
 
        
        for k = 1:L-1
            observ(j,datal(k)+1) = observ(j,datal(k)+1) + 1;
        end
        Nq = max(Nq,max(datal(1:length(datal)-1)));
    end
    fclose(f1);
    observ = observ(:,1:Nq+1);
    V = observ;
    m = size(observ,1);
    for i = 1:m
        for j = 1:Nq
            observ(i,1) = observ(i,1) * (1/sqrt(2*w(j)))^(observ(i,j+1));
        end
    end
    
    for k = 1:Nq
        m = m + 1;
        observ(m,1) = 1/2 * w(k) / 2;
        observ(m,k+1) = -2;
        V(m,1) = 1/2;
        V(m,k+1) = -2;
    end
    
    
end
