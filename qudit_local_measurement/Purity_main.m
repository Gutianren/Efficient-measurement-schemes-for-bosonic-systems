global Nq meas pr len observ d phi expect coeff purity_expect rho
d = 4;
Nq = 3;
file = 'H2OCoe0';
[expect, phi] = GetGroundEnergy(d,file);
p = 0.05;
rho_id = eye(d^Nq) / d^Nq;
rho = (1-p) * phi * phi' + p * rho_id;
purity_expect = Purity(rho);

Purity_Getobserv();
Purity_Getcoefficient();
Purity_OGM();
Purity_Sample_main();

