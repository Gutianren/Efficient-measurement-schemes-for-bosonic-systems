function main(dd,TT,Repeat)
    %clear all;
    global d file expect phi pr Nq observ meas len
    d = dd;
    n = 3;
    file = 'H2OCoe0';
    %[expect,phi] = GetGroundEnergy(d,file);
    %H = GetGroundEnergy(d,file);
    [H,w] = GetHamiltonian(d,file);
    [expect,phi] = AnyPhi(d,n,H);
    tic;
    diag_var = CutOGM(d,file,0);
    tOGM = toc;
    m1 = size(observ,1);
    l1 = size(meas,2);
    tic;
    var = varianceTotal(d);
    tvar = toc;
    fprintf('\nvar = %f\n', var);

    for t = 1:length(TT)
        [expect,phi] = AnyPhi(d,n,H);
        T = TT(t);
        tic;
        [ex,err,data] = Sample_main(0,T,Repeat);
        tsample = toc;
        
        format = './zdata_d_ghz_d%dR%dT%d.txt';
        outfile = sprintf(format,d,Repeat,T);
        save(outfile,'m1','data','ex','expect','err','-ascii');
    end
end

% phi = BinaryGroundState(d,phi);
% 
% tic;
% diag_var_qubit = CutOGM(d,file,1);
% tOGM_qubit = toc;
% m2 = size(observ,1);
% l2 = size(meas,2);
% tic;
% var_qubit = varianceTotal(2);
% tvar_qubit = toc;
% fprintf('\nbitvar = %f\n',var_qubit);
% tic;
% [ex_qubit, err_qubit] = Sample_main(1);
% tsample_qubit = toc;    
% GetObservableExpect()
 
% fprintf('------------Final Result for d = %d---------------\n',d);
% fprintf('Exact Expectation = %f\n', expect);
% fprintf('---Qudit---\n');
% fprintf('Encoded observables: %d, Sets: %d\n', m1, l1);
% fprintf('Diagonal variance: %f, Total variance: %f\n', diag_var, var);
% fprintf('Sample Expectation: %f, Sample Error: %f\n', ex, err);
% fprintf('Time: (OGM) %f, (varicance) %f, (sample) %f\n', tOGM, tvar, tsample);
% fprintf('---Qubit---\n');
% fprintf('Encoded observables: %d, Sets: %d\n', m2, l2);
% fprintf('Diagonal variance: %f, Total variance: %f\n', diag_var_qubit, var_qubit);
% fprintf('Sample Expectation: %f, Sample Error: %f\n', ex_qubit, err_qubit);
% fprintf('Time: (OGM) %f, (varicance) %f, (sample) %f\n', tOGM_qubit, tvar_qubit, tsample_qubit);
