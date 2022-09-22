function [coes, strs] = transform_compact(d,mat)
    % transform a d*d matrix into len = ceil(log2(d)) Pauli strings
    K = ceil(log2(d));
    coes = [];
    strs = [];
    len = 0;
    coeV = [1,1;1,-1i;1,1i;1,-1];
    strV = [0,3;1,2;1,2;0,3];
    for i = 1:d
        for j = 1:d
            if mat(i,j) ~= 0
                coe = [mat(i,j)/2^K];
                str = []; 
                bis = dec2bin(i-1,K);
                bjs = dec2bin(j-1,K);
                bi = str2num(bis(:));
                bj = str2num(bjs(:));
                for k = 1:K
                    x = bi(k) + 2*bj(k) + 1;
                    coev = coeV(x,:)';
                    strv = strV(x,:)';
                    [coe,str] = strmul(coe,str,coev,strv);
                end
                [coes, strs] = mergestr(coes,strs,coe,str);
            end
        end
    end
end