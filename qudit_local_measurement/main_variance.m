global d file expect phi
d = 3;
file = 'H2OCoe0';
[expect,phi] = GetGroundEnergy(d,file);
diag_var = CutOGM(d,file,0);
var = varianceTotal(d);
fprintf('\nvar = %f\n', var);

% next we deal with the qubit version
phi = BinaryGroundState(d,phi);
diag_var_bit = CutOGM(d,file,1);
var = varianceTotal(2);
fprintf('\nbitvar = %f\n',var);


