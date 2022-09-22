function [coes, strs] = encoding(d,file)
    % This function encodes given molecular potential into d-dimension GGB
    % representations. coes are coefficients and strs are n 1~d*d numbers
    % denoting the matrix

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
        
        [coeH, strH] = decomposestr(d,n,H);
        coeH = coeH * k;
        [coes,strs] = mergestr(coes,strs,coeH,strH);
        %coes = [coes;coeH];
        %strs = [strs;strH];
    end
    fclose(f);
    
    % encoding the kinetic energy p^2/2
    H = zeros(d,d,n);
    for i = 1 : n
        H(:,:,i) = eye(d);
    end
    for i = 1:n
        H(:,:,i) = - p0 * p0 * w(i)/2 / 2;
        [coeH, strH] = decomposestr(d,n,H);
        % coeH = coeH * k;
        [coes,strs] = mergestr(coes,strs,coeH,strH);
        
%         coes = [coes;coeH];
%         strs = [strs;strH];
        H(:,:,i) = eye(d);
    end
    
%     % merge the same terms
%     newcoes = [];
%     newstrs = [];
%     for i = 1:length(coes)
%         f = 0;
%         for j = 1:length(newcoes)
%             if newstrs(j,:) == strs(i,:)
%                 newcoes(j) = newcoes(j) + coes(i);
%                 f = 1;
%                 break;
%             end    
%         end
%         if f == 0
%             newcoes = [newcoes;coes(i)];
%             newstrs = [newstrs;strs(i,:)];
%         end
%     end
%     coes = newcoes;
%     strs = newstrs;
    
    % change Identity to zero
    for i = 1:length(coes)
        for j = 1:n
            if strs(i,j) == d*d
                strs(i,j) = 0;
            end
        end
    end
    
    fprintf('Total %d observables encoded\n', length(coes));
end