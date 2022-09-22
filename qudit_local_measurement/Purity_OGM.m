function Purity_OGM()
    
    % encoding
    global Nq meas pr len observ d coeff
    
    MAX_step = 10;
 
    [coeff, observ] = sort_observ(coeff,observ); %sorting
    observ(:,1) = coeff;
    
    m = size(observ,1);
    
    % Grouping
    cur = 0; % Number of sets
    added = zeros(1,m);
    
    meas = [];
    while 1
        j = 1;
        while j <= m && added(j) > 0
           j = j + 1; 
        end
        if j > m 
           break; 
        end
        cur = cur + 1;
        ele_cur(cur) = 1; % current element index in the current set
        meas(:,cur)=observ(j,2:Nq+1); % P initialized to Q
        set(ele_cur(cur),cur)=j; % set(i,j)=k, the i-th elements in set j is the k-th elements in observ.
        added(j) = 1;
        
        sum_pr = abs(observ(j,1));%initialize pr
        start = j;
        j = j + 1; %% current observable
        
        %% initialize the set elements 
        while j <= m 
            if IfCommute(meas(:,cur),observ(j,2:Nq+1))
                ele_cur(cur) = ele_cur(cur) + 1;
                set(ele_cur(cur), cur) = j;
                meas(:,cur) = update_meas(meas(:,cur),observ(j,2:Nq+1));
                added(j) = 1;
                sum_pr = sum_pr + (abs(observ(j,1)));
            end
            j = j + 1;
        end
        
        %% initize pr
        pr(cur) = sum_pr;
        %%maximize the set
        for k = 1 : start - 1
            if IfCommute(meas(:,cur),observ(k,2:Nq+1))
                ele_cur(cur) = ele_cur(cur) + 1;
                set(ele_cur(cur),cur) = k;
                meas(:,cur) = update_meas(meas(:,cur), observ(k,2:Nq+1));
            end
        end      
    end
    
    %% the number of sets.
    len = cur;

    format = 'The number of sets: %d.\n';
    fprintf(format, len);
    
    pr = pr / sum(pr);
    
    tic
    exitflag = 0;
    diag_new = 0;
    T = 1000; %% if the result is not good, we should amplify this.
    for step = 1 : MAX_step
        %%sorting the meas with probability descending.
        measAndP(1:Nq, 1:len) = meas(:,1:len);
        measAndP(Nq+1, 1:len) = pr(1:len);
        measAndP = measAndP(:,1:len);
        measAndP = sort_meas(measAndP);%sort
        meas = measAndP(1:Nq,1:len);
        pr = measAndP(Nq+1,1:len);
    % % % %     [new_len, new_pr] = CutDownSet(pr, T);%%cutdown
    % % % %     if abs(new_len - len) < Nq
    % % % %        break; 
    % % % %     end
       [new_len, new_pr] = CutMoreSet(pr, T);%%cutmore
       if variance(new_pr) < variance(pr)
           pr = new_pr;
           len = new_len;
           meas = meas(:, 1:len);
       elseif exitflag == 1
          break; 
       end
        [pr, diagonal_var, exitflag] = OptDiagVar(pr, 10 * step);
        if abs(diag_new - diagonal_var)<1e-3
           break; 
        end
        diag_new = diagonal_var;
       format = '%d-th round: The number of sets: %d, diag_variance: %f.\n';
       fprintf(format, step, len, diagonal_var);
        
    end
    if exitflag == 0
        [pr, diagonal_var, exitflag] = OptDiagVar(pr, 100);
    end

    timeSpend = toc;
    fprintf('CutOGM, time cost:%f; %d rounds in total; Optimal diagonal var: %f\n', timeSpend, step, diagonal_var);

end
