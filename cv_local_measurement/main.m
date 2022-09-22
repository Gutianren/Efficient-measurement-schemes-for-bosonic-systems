function main(dd,T,Repeat)
    %clear all;
    global Nq observ e gs w m meas pr d phi phi_p V x p len expect dx dp
    file = 'H2OCoe0';
    d = dd;
    %T = 1000;
    %Repeat = 100;
    GetObserv(file);
    CutOGM();
    %file1 = ['gs_',file,'_', num2str(d), '.mat'];
    %load(file1);
    [H,w] = GetHamiltonian(d,'H2OCoe0');
    [e,gs] = AnyPhi(d,Nq,H);
    m1 = 2;
    dx = 1;
    xmax = 150;
    
    x = -xmax:dx:xmax;
    
    dp = 0.01;
    pmax = 1;
    p = -pmax:dp:pmax;
    
    [phi,phi_p] = GetGroundState(x,p);
    normx = sum(abs(phi).^2)*dx;
    normp = sum(abs(phi_p).^2)*dp;
    
    if max(abs(normx-1)) < 1e-10 && max(abs(normp-1)) < 1e-10
        expect = GetExpectation(x,p)
        [ex,err,data] = Sample_main(T,Repeat);
    else
        fprintf('Improper x,p range\n');
        error_x = abs(normx-1)
        error_p = abs(normp-1)
    end
        
    format = './zdata_c_ghz_d%dR%dT%d.txt';
    outfile = sprintf(format,d,Repeat,T);
    save(outfile,'m1','data','ex','expect','err','-ascii');
end