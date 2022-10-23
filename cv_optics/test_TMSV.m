
global rho_x rho_V T Repeat ex error alpha Vidx Nq

Nq = 2;
alpha = [-1/2,1/4,1/4];
Vidx = [0,0;2,0;-2,0]

Repeat = 100;
r_ex = [];
r_err = [];
rr = 0:0.05:2;
for r = rr
    T = 1000;
    [rho_x,rho_V] = TMSV(r);
    main();
    r_ex = [r_ex,ex];
    r_err = [r_err,error];
end
r_save = [rr', r_ex', r_err'];
outfile = 'TMSV_r.txt';
save(outfile,'r_save','-ascii');

Repeat = 100;
r = 1;
T_ex = [];
T_err = [];
TT = 200:200:10000;
for T = TT
    [rho_x,rho_V] = TMSV(r);
    main();
    T_ex = [T_ex,ex];
    T_err = [T_err,error];
end
T_save = [TT', T_ex', T_err'];
outfile = 'TMSV_T.txt';
save(outfile,'T_save','-ascii');