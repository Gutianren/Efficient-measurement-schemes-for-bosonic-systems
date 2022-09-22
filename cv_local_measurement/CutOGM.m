function CutOGM()
    global observ len pr Nq m meas len
    
    observ = sort_observ(observ);
    
    cur = 0;
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
        %ele_cur(cur) = 1; % current element index in the current set
        meas(:,cur) = sign(observ(j,2:Nq+1)); % P initialized to Q
        %set(ele_cur(cur),cur)=j; % set(i,j)=k, the i-th elements in set j is the k-th elements in observ.
        added(j) = 1;
        
        sum_pr = abs(observ(j,1));%initialize pr
        start = j;
        j = j + 1; %% current observable
        
        %% initialize the set elements 
        while j <= m 
            if IfCommute(meas(:,cur),observ(j,2:Nq+1))
                %ele_cur(cur) = ele_cur(cur) + 1;
                %set(ele_cur(cur), cur) = j;
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
                %ele_cur(cur) = ele_cur(cur) + 1;
                %set(ele_cur(cur),cur) = k;
                meas(:,cur) = update_meas(meas(:,cur), observ(k,2:Nq+1));
            end
        end      
    end
    
    %% the number of sets.
    len = cur;
    pr = pr / sum(pr);
    
    
end
