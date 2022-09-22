function [trace] = TraceGlobal(d,phi,Q1,Q2)
    % Calculate trace(|phi><phi|*Q1*Q2)   
    % Here Q1,Q2 can be string of GGB matrixes
    % d is the dimension of each qubit
    % phi is a state of d^n coefficients, denoting |0..0> to |d-1,..,d-1>
    
    dimension = length(phi);
    Nq = floor(log(dimension)/log(d));
    
    if Nq ~= length(Q1) || Nq ~= length(Q2)
        disp('The qubits of the input state does not match the dimension of Q.\n');
        return;
    end
    
    coeff = 1; % to extract the overall coefficient
    reQ1 = zeros(1,Nq); % to extract the plain
    reQ2 = zeros(1,Nq);
    for qubit = 1:Nq
        if abs(Q1(qubit)) > 0
            reQ1(qubit) = floor(abs(Q1(qubit)));
            coeff = coeff * Q1(qubit) / reQ1(qubit);
        end
        if abs(Q2(qubit)) > 0
            reQ2(qubit) = floor(abs(Q2(qubit)));
            coeff = coeff * Q2(qubit) / reQ2(qubit);
        end
    end

    k = 0;
    for j = 1 : dimension
        if phi(j) ~= 0
            k = k + 1;
            s(k) = j; %%s: contains all of elements with non-zero values.
        end
    end
    
%     if sum(reQ==0) < Nq % containing non-diagonal terms
        trace = 0;
        for j = 1:k
            target = CalculateTarget(d,s(j),Q1,Q2); % Act Q on s(j)
            if ~isempty(find(s==target, 1))
                trace = trace + real(phi(s(j))' * phi(target) * TraceLocalGen(d, s(j),target, Q1, Q2));
            end
        end
%       
%     else % only diagonal terms
%         trace = 0;
%         for j = 1 : k
%             trace = trace + phi(s(j))' * phi(s(j)) * TraceLocalGen(s(j),s(j),Q1, Q2); %<j|Q|l>
%         end        
%     end
    
end
    