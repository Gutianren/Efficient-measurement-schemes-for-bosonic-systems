function [coes,strs]=decomposestr_qubit(d,n,H,mode)
    % For given local Hamiltonian(in d*d*n tensor), decompose it to
    % coefficient-strings
    % mode == 0 : Direct mapping, mode == 1: Compact mapping
    
    if mode == 0
        [coes,strs] = transform(d,H(:,:,1));
        for i = 2:n
            [coe,str] = transform(d,H(:,:,i));
            [coes,strs] = strmul(coes,strs,coe,str);
        end
    else
        [coes,strs] = transform_compact(d,H(:,:,1));
        for i = 2:n
            [coe,str] = transform_compact(d,H(:,:,i));
            [coes,strs] = strmul(coes,strs,coe,str);
        end
    end
    
end