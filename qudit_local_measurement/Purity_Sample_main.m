function [ex, error] = Purity_Sample_main()

    global observ meas len phi pr rho purity_expect

    T = 1000;
    Repeat = 100;
    
    [m, Nq] = size(observ);
    Nq = Nq - 1;
    format = 'The number of qubits: %d, the number of observables: %d';
    fprintf(format, Nq, m);

    %%% convert phi to s, where Coord_phi_t = j represents the j-th coordinate equals 1.
%     k  = 0;
%     dimension = length(phi);
%     for j = 1 : dimension
%         if phi(j) ~= 0
%             k = k + 1;
%             Coord_phi(k,1) = j-1; %%s: contains all of elements with non-zero values, start from 0-->0...0
%             Coord_phi(k,2) = phi(j);
%         end
%     end
%     % fprintf('The number of non-zero basis: %d\n', k);
    %%Sorting the meas and pr.

    measAndP(1:Nq, 1:len) = meas(:,1:len);
    measAndP(Nq+1, 1:len) = pr(1:len);
    measAndP = measAndP(:,1:len);

    measAndP = sort_meas(measAndP);
    len = size(measAndP,2);
    meas = measAndP(1:Nq,:);
    pr = measAndP(Nq + 1,:);
    j = 0;
    s = 0;
    
    while j < T
        if s == len
            break;
        end
        s = s + 1;
        % NumMeas = ceil(T * pr(s));
        rx = rand();
        amount = T * pr(s);
        if amount > 1
            NumMeas = floor(amount);
            if rx < mod(amount, NumMeas)
                NumMeas = NumMeas + 1;
            end
        else
            NumMeas = 1;
        end
        for k = 1 : NumMeas
            j = j + 1;
            measure(j,:) = meas(:,s).';
            if j>= T
                break;
            end
        end
    end
    len = s;
    meas = meas(:,1:len);
    pr = pr(1:len);
    pr = pr / sum(pr);

    %%upadate pr
%     sum = 0;
%     for j = 1 : len
%         sum = sum + pr(j);
%     end
%     for j = 1 : len
%         pr(j) = pr(j)/sum;
%     end

    T = size(measure,1);
    fprintf('The number of sampling: %d\n', T);
    sum_err = 0;
    ex = 0;
    tic;
    for k = 1 : Repeat
        v = 2 * Purity_Sample_DerandOGM(observ, Nq, m, measure, T, rho) - 1;
        ex = ex + v;
        sum_err = sum_err + (v-purity_expect)^2;
        fprintf('Iteration %d of %d, v=%f\n',k,Repeat,v);
    end
    ex = ex / Repeat;
    error = sqrt(sum_err/Repeat);
    timeElapsed = toc;

    fprintf('expect_purity = %f, error = %f\n', ex, error);
    fprintf('Expectation is %f, total run time: %f\n', purity_expect, timeElapsed);
end