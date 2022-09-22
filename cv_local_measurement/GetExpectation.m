function [ex] = GetExpectation(x,p)
    global phi phi_p gs d Nq V w
    ex = 0;
    dx = x(2)-x(1);
    dp = p(2)-p(1);
    L = length(x);
    
    for l = 1:size(V,1) 
        for j1 = 1:d^Nq
            for j2 = 1:d^Nq
                v1 = ToDinary(d,j1-1,Nq); % The state is gs(j) * phi1_v1[x1] * phi2_v2[x2] * phi3_v3[x3]
                v2 = ToDinary(d,j2-1,Nq);
                prod = gs(j1)' * gs(j2) * V(l,1);
                for k = 1:Nq % k-th particle
                    if V(l,k+1) > 0 % Potential
                        prod = prod * (x.^V(l,k+1)) * (conj(phi(:, (k-1)*d+v1(k)+1)) .* phi(:, (k-1)*d+v2(k)+1)) * dx;
                    else % Kinetic
                        prod = prod * (p.^(-V(l,k+1))) * (conj(phi_p(:, (k-1)*d+v1(k)+1)) .* phi_p(:, (k-1)*d+v2(k)+1)) * dp;                    
                    end
                end
                ex = ex + prod;
            end
        end
    end
    
%     ex1 = 0;
%     for j = 1:d^Nq
%         v = ToDinary(d,j-1,Nq);
%         pd1 = 0;
%         for k = 1:Nq
%             pd1 = pd1 + w(k) * (v(k) + 1/2);
%         end
%         ex1 = ex1 + pd1 * abs(gs(j))^2;
%     end
%     fprintf('ex1=%f',ex1);
    
    % Kinetic energy
%     for l = 1:Nq
%         for j = 1:d^Nq
%             v = ToDinary(d,j-1,Nq);
%             prod = gs(j)^2 * 1/2;
%             for k = 1:Nq
%                 if l == k
%                     prod = prod * (p.^2) * (abs(phi_p(:, (k-1)*d+v(k)+1)).^2) * dp;
%                 end
%             end
%             ex = ex + prod;
%         end
%     end        
%         
end