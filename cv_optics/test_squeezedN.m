
global rho_x rho_V T Repeat ex error alpha Vidx Nq mode expect


Repeat = 100;
M = 500;
%N_state = 10;
T0 = 1000;

%Nq = 10;
K = [2,4,6];

result = [];
for Nq = 2:16
    beta = (ones(Nq)-eye(Nq))/Nq;    % equally squeezed
    rho_V = squeezedN(beta);    
    rho_x = zeros(2*Nq,1);  % without displacement
    
    for k = K
        for j = 1:5
            alpha = rand(1,M);
            alpha = alpha / norm(alpha);        
            Vidx = [];
            for i = 1:M
                Vidx(i,:) = rndpx(Nq,k);
            end
        
            for mode = 3:3
                fprintf('Currently Running Nq=%d, k=%d, j=%d, mode=%d',Nq,k,j,mode)
                T = T0;
                main();
                result = [result;[Nq,k,mode,expect,ex,error]];
            end
        end
    end
end
result = real(result);
yogm = result(find(result(:,3)==1),6);
yldf = result(find(result(:,3)==2),6);
yis = result(find(result(:,3)==3),6);
ycs = result(find(result(:,3)==4),6);
x = result(find(result(:,3)==4),1);
save RAN_k_squeezed_n2:16_onlyIS.mat