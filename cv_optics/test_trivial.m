

global rho_x rho_V T Repeat ex error alpha Vidx Nq mode

mode = 1;

Nq = 2;
alpha = [-1/2,1/4,1/4];
Vidx = [0,0;2,0;-2,0]

Repeat = 100;
T = 1000;
r=1;
[rho_x,rho_V] = TMSV(r);
main();
