function [ex,error,results] = Sample_main(T,Repeat)
    
    global observ meas phi phi_p x p m Nq pr len expect
    
    %T = 1000;
    %Repeat = 100;
    format = 'The number of qubits: %d, the number of observables: %d\n';
    fprintf(format, Nq, m);

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
    
    tic;
    T = size(measure,1);
    results = zeros(Repeat,2);
    fprintf('The number of sampling: %d\n', T);
    sum_err = 0;
    ex = 0;
    for k = 1:Repeat
        v = Sample_DerandOGM(measure, T);
        ex = ex + v;
        sum_err = sum_err + (v-expect)^2;
        fprintf('Iteration %d of %d, v=%f\n',k,Repeat,v);
        results(k,:) = [v,sqrt((v-expect)^2)];   
    end
    ex = ex / Repeat;
    error = sqrt(sum_err/Repeat);
    timeElapsed = toc;

    fprintf('expect_m = %f, error = %f\n', ex, error);
    fprintf('Expectation is %f, total run time: %f\n',expect, timeElapsed);
    

end