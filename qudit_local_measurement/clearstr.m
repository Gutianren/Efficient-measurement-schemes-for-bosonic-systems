function [newcoes,newstrs] = clearstr(coes,strs)
    % merge the same terms
    
    ncoes = [];
    nstrs = [];
    for i = 1:length(coes)
        f = 0;
        for j = 1:length(ncoes)
            if nstrs(j,:) == strs(i,:)
                ncoes(j) = ncoes(j) + coes(i);
                f = 1;
                break;
            end    
        end
        if f == 0
            ncoes = [ncoes;coes(i)];
            nstrs = [nstrs;strs(i,:)];
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