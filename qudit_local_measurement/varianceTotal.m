function var = varianceTotal(d)
    % d : qudit
    % d == 2: qubit
    
    global observ meas len phi expect pr
    
    %expect
    var = 0;
    m = size(observ,1);
    n = size(observ,2) - 1;
    
    for j = 1 : m
        for k = 1 : m
            temp_j = 0;
            temp_k = 0;
            temp_n = 0; % numerator
            temp_d = 0; % denominator
            for i = 1:len
                flag_1 = 0; 
                flag_2 = 0;
                if IfCommute(observ(j, 2:n+1), meas(:,i))
                    temp_j = temp_j + pr(i);
                    flag_1 = 1;
                end
                if IfCommute(observ(k, 2:n+1), meas(:,i))
                    temp_k = temp_k + pr(i);
                    flag_2 = 1;
                end
                if flag_1 == 1 && flag_2 == 1
                    temp_n = temp_n + pr(i);
                end    
            end
            temp_d = temp_j * temp_k;
            % g_factor = temp_n / temp_d;
            % [~,mul] = strmul(1, observ(j,2:n+1), 1, observ(k,2:n+1));
            % since GGB matrixes do not form a group
            % We need the trace here(need the ground state)
            if temp_d ~= 0
                trace = TraceGlobal(d,phi,observ(j,2:n+1),observ(k,2:n+1));
                var = var + observ(j,1) * observ(k,1) * temp_n/temp_d * trace;
            end
        end
        fprintf('\nRound %d of %d', j,m);
    end
    
    %diag_var = variance(pr);
    %var = 2 * var +  diag_var;
    var = var - expect ^ 2;

end
