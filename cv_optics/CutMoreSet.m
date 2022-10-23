function [len, pr] = CutMoreSet(pr, T)
%% T *p1, ..., T*ps, cut measurements with indices: s+1..length(pr)


    count = 0;
    s = 0;
    len = length(pr);

    while count < T && s < len
        s = s + 1;
        %%no rand here
        count = count + ceil(T * pr(s));
    end

    if s < len
        len = s;

        pr = pr(1:len);
        
        pr = pr / sum(pr);
    end
end