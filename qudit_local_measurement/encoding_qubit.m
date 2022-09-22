function [coes, strs] = encoding_qubit(d, file, mode)
    % encoding to qubit results
    % mode == 0 : non-compact, mode == 1 : compact
    % For the lowest d energies
    
    [ad,a] = geta(d);
    q0 = ad + a;
    p0 = ad - a;
    
    n = 0;
    file = strcat(file,'.txt');
    f1 = fopen(file);
    while ~feof(f1)
        fline = fgetl(f1);
        datal = str2num(fline);
        n = max(n,max(datal(1:length(datal)-1)));
    end
    fclose(f1);
    
    coes = [];
    strs = [];
    % encoding potentials
    f = fopen(file);
    while ~feof(f)
        fline = fgetl(f);
        datal = str2num(fline);
        j = length(datal)-1;
        idx = datal(1:j);
        k = datal(j+1);
        
        if j == 2 && idx(1) == idx(2)
           i = idx(1);
           w(i) = sqrt(2*k);
        end
        
        H = zeros(d,d,n);
        for i = 1 : n
            H(:,:,i) = eye(d);
        end
        for i = 1 : j
            H(:,:,idx(i)) = H(:,:,idx(i)) * q0 / sqrt(2*w(idx(i)));
        end
        
        [coeH, strH] = decomposestr_qubit(d,n,H,mode);
        coeH = coeH * k;
        %coes = [coes;coeH];
        %strs = [strs;strH];
        [coes,strs] = mergestr(coes,strs,coeH,strH);
    end
    fclose(f);
    
    % encoding the kinetic energy p^2/2
    H = zeros(d,d,n);
    for i = 1 : n
        H(:,:,i) = eye(d);
    end
    for i = 1:n
        H(:,:,i) = - p0 * p0 * w(i)/2 / 2;
        [coeH, strH] = decomposestr_qubit(d,n,H,mode);
        %coeH = coeH * k;
        %coes = [coes;coeH];
        %strs = [strs;strH];
        [coes,strs] = mergestr(coes,strs,coeH,strH);
        H(:,:,i) = eye(d);
    end
    
    %[coes,strs] = clearstr(coes,strs);
    fprintf('Total %d observables encoded\n', length(coes));
end
    
        