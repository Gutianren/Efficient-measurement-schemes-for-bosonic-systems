function CS()
    global observ len pr Nq m meas T
    
    pr = ones(1,2^Nq);
    pr = pr / sum(pr);
    
    meas = zeros(Nq,2^Nq);
    meas(:,1) = ones(Nq,1);
    for j = 2:2^Nq
        meas(:,j) = meas(:,j-1);
        k = Nq;
        while meas(k,j) == -1
            meas(k,j) = 1;
            k = k - 1;
        end
        meas(k,j) = -1;
    end
    len = 2^Nq;
end
