
global Nq m len rho_x rho_V observ expect pr meas mode

%Nq = 2;
%M = 2*Nq;
%len = 2;
%alpha = 1:2*Nq;
%Vidx = kron(eye(Nq),[1;-1]);


%alpha = 1;
%Vidx = zeros(1,Nq);
%Vidx(1) = 4;


observ = [alpha',Vidx];

%observ = zeros(1,21);
%observ(1) = 1;
%observ(2) = 2;
%observ(3) = -2;

m = size(observ,1);
%rho_x = rand(1,2*Nq)';

%rho_x = (1:2*Nq)';
%rho_x = zeros(2*Nq,1);

%A = rand(2*Nq);
%rho_V = A * A' / (2*Nq);

%rho_V = eye(2*Nq); % Thermal State

%[rho_x,rho_V] = TMSV(1);

expect = GetExpectation(rho_x',rho_V,observ)

if mode == 1
    CutOGM();
elseif mode == 2
    LDFGroup();
elseif mode == 3
    ImportanceSampling();
elseif mode == 4
    CS();
end

%% Simple Important-Sampling Process

% meas = [ones(Nq,1),-ones(Nq,1)];
% px = sum(alpha(1:Nq));
% pp = sum(alpha(Nq+1:2*Nq));
% norm = abs(px) + abs(pp);
% px = px / norm;
% pp = pp / norm;
% pr = [px,pp];
% pr = [0.5,0.5]

%%

Sample_main();

