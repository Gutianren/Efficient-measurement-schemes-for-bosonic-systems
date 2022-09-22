global d file expect phi pr Nq observ meas len
d = 4;
file = 'H2OCoe0';
[expect,phi] = GetGroundEnergy(d,file);
phi = BinaryGroundState(d,phi);
diag_var = CutOGM(d,file,1);
Sample_main(1);
