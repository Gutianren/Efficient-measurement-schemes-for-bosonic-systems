function [ex, error, results] = Sample_main(mode,T,Repeat)

    global observ meas len phi pr expect

    %T = 1000;
    %Repeat = 100;
    
    [m, Nq] = size(observ);
    Nq = Nq - 1;
    format = 'The number of qubits: %d, the number of observables: %d';
    fprintf(format, Nq, m);

    %prompt = 'Please input a number between 0 and 6 to represent the program you wish to execute.\n(0: OGM; 1: LBCS; 2: grouping; 3: Derandomized alg for Huang et al; 4: Derandomized alg for OGM;):\n';
%     str = 'Hamiltonian/';
% 
%     mole = file;
%     file = strcat(file,'.txt');
%     strF = strcat(str, file);
% 
% 
%     observ = load(strF);
%     [m,Nq] = size(observ);
%     Nq = Nq - 1;
%     format = 'The number of photons: %d, the number of observables: %d';
%     fprintf(format, Nq, m);
%     display(strF);
% 
%     phi = get_qstateGlobal(mole); %%mole = file (e.g.: LiH_12jw)
% 
% 
%     %expect = energy(phi, observ); replace computing to read from file
%     str_energy = 'energy/ground_state_';
%     str_energy = strcat(str_energy, mole);
%     str_energy = strcat(str_energy, '.txt');
%     f_energy = fopen(str_energy,'r');
%     temp = fgetl(f_energy);
%     [energy_mole]= sscanf(temp,'%s %f');
%     %expect = -28.3565
%     long_energy= length(energy_mole);
%     expect = energy_mole(long_energy);
% 
%     tic
% 
%     if x == 1
%         str_OGM = 'CutSet/OGMV1_';
%     elseif x == 2
%         str_OGM = 'CutSet/OGMV2_';
%     else
%         disp('Error');
%         return;
%     end
% 
%     str_OGM = strcat(str_OGM, file);
%     measAndP = load(str_OGM);
    
    %%% convert phi to s, where Coord_phi_t = j represents the j-th coordinate equals 1.
    k  = 0;
    dimension = length(phi);
    for j = 1 : dimension
        if phi(j) ~= 0
            k = k + 1;
            Coord_phi(k,1) = j-1; %%s: contains all of elements with non-zero values, start from 0-->0...0
            Coord_phi(k,2) = phi(j);
        end
    end
    fprintf('The number of non-zero basis: %d\n', k);
    %%Sorting the meas and pr.

    measAndP(1:Nq, 1:len) = meas(:,1:len);
    measAndP(Nq+1, 1:len) = pr(1:len);
    measAndP = measAndP(:,1:len);

    measAndP = sort_meas(measAndP);
    len = size(measAndP,2);
    meas = measAndP(1:Nq,:);
    pr = measAndP(Nq + 1,:);
    j = 0;
    s = 0;
    
    while j < T
        if s == len
            break;
        end
        s = s + 1;
        % NumMeas = ceil(T * pr(s));
        rx = rand();
        amount = T * pr(s);
        if amount > 1
            NumMeas = floor(amount);
            if rx < mod(amount, NumMeas)
                NumMeas = NumMeas + 1;
            end
        else
            NumMeas = 1;
        end
        for k = 1 : NumMeas
            j = j + 1;
            measure(j,:) = meas(:,s).';
            if j>= T
                break;
            end
        end
    end
    len = s;
    meas = meas(:,1:len);
    pr = pr(1:len);
    %%upadate pr
    sum = 0;
    for j = 1 : len
        sum = sum + pr(j);
    end
    for j = 1 : len
        pr(j) = pr(j)/sum;
    end
    
    tic;
    results = zeros(Repeat,2);
    T = size(measure,1);
    fprintf('The number of sampling: %d\n', T);
    sum_err = 0;
    ex = 0;
    for k = 1 : Repeat
        if mode == 0
            v = Sample_DerandOGM(observ, Nq, m, measure, T, Coord_phi);
        else
            v = Sample_DerandOGM_qubit(observ, Nq, m, measure, T, Coord_phi);
        end
        ex = ex + v;
        results(k,:) = [v,sqrt((v-expect)^2)];
        sum_err = sum_err + (v - expect)^2;
        fprintf('Iteration %d of %d, v=%f\n',k,Repeat,v);
    end
    ex = ex / Repeat;
    error = sqrt(sum_err/Repeat);
    timeElapsed = toc;

    fprintf('expect_m = %f, error = %f\n', ex, error);
    fprintf('Expectation is %f, total run time: %f\n',expect, timeElapsed);
end