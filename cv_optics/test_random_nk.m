
global rho_x rho_V T Repeat ex error alpha Vidx Nq mode

%Nq = 10;
Repeat = 100;
M = 100;
N_state = 10;
T0 = 1000;

NQ = 13:20;
ran_nk = [];


for mode = 1:4

    for Nq = NQ
        
        clear RV;
        for i = 1:N_state
            A = rand(2*Nq);
            RV(:,:,i) = A * A';
            RV(:,:,i) = RV(:,:,i) / trace(RV(:,:,i));
        end
        
        
        %kmax = 0;
        K = [2,4,6];
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
                ran_nk = [ran_nk;[Nq,kmax,ex,error]];
            end
        end
    end
end
save RAN_nk_3.mat;
outfile = 'RAN_nk_3.txt';
save(outfile,'ran_nk','-ascii');
