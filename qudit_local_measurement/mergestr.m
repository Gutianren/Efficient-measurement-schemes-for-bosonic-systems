function [newcoes,newstrs] = mergestr(coes,strs,coes1,strs1)
    % merge the same terms
    ncoes = coes;
    nstrs = strs;
    for i = 1:length(coes1)
        f = 0;
        for j = 1:length(coes)
            if nstrs(j,:) == strs1(i,:)
                ncoes(j) = ncoes(j) + coes1(i);
                f = 1;
                break;
            end    
        end
        if f == 0
            ncoes = [ncoes;coes1(i)];
            nstrs = [nstrs;strs1(i,:)];
        end
    end
    newcoes = [];
    newstrs = [];
    for i = 1:length(ncoes)
        if ncoes(i) ~= 0
            newcoes = [newcoes;ncoes(i)];
            newstrs = [newstrs;nstrs(i,:)];
        end
    end
    
end