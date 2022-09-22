function [coes, strs] = decomposestr(d,n,mats)
    % decompose a Hermitian string to M GGB strings
    % each string has n elements
    [coes, strs] = decompose(d,mats(:,:,1));
    M = length(coes);
    for i = 2:n
        [coe, str] = decompose(d,mats(:,:,i));
        [coes,strs] = strmul(coes,strs,coe,str);
    end
end