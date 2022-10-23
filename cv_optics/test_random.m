
global rho_x rho_V T Repeat ex error alpha Vidx Nq

Nq = 10;
Repeat = 100;
M = 100;
N_state = 10;
T0 = 1000;


for i = 1:N_state
    %RX(:,i) = rand(2*Nq,1);
    %RX(:,i) = RX(:,i) - (sum(RX(:,i)) / Nq);
    A = rand(2*Nq);
    RV(:,:,i) = A * A';
    RV(:,:,i) = RV(:,:,i) / trace(RV(:,:,i));
end

ran_k = [];

%kmax = 0;
K = 2:12;
for kmax = K
    
    alpha = rand(1,M);
    alpha = alpha / norm(alpha);

    Vidx = [];
    for i = 1:M
        kk = kmax;
        %kk = randi(kmax);
        Vidx(i,:) = rndpx(Nq,kk);
    end
    ran_ex = [];
    ran_err = [];
    for i = 1:N_state
        %rho_x = RX(:,i);
        rho_x = zeros(2*Nq,1);
        rho_V = RV(:,:,i);
        T = T0;
        main();
        ran_k = [ran_k;[kmax,ex,error]];
    end
end
outfile = 'RAN_k_1.txt';
save(outfile,'ran_k','-ascii');
