%% Sensor Selection Example
warning off
%% Clear
clear, close all
SWA = 'uswa';  % running with uncertainty and noise
% SWA = 'swa';  % running without uncertainty and noise
addpath('../../lib/')

%% Model Parameters for 5 Scenarios
model(1).C = eye(6);
model(1).D = [0;0;0;0;0;0];
model(2).C = [0 1 0 0 0 0;
              0 0 1 0 0 0;
              0 0 0 1 0 0;
              0 0 0 0 1 0;
              0 0 0 0 0 1];
model(2).D = [0;0;0;0;0];
model(3).C = [0 0 1 0 0 0;
              0 0 0 1 0 0;
              0 0 0 0 1 0;
              0 0 0 0 0 1]; % This row is added because of a bug in Yalmip
model(3).D = [0;0;0;0];
model(4).C = [1 0 0 0 0 0;
              0 1 0 0 0 0;
              0 0 1 0 0 0;
              0 0 0 1 0 0]; % This row is added because of a bug in Yalmip
model(4).D = [0;0;0;0];
model(5).C = [0 1 0 0 0 0;
              0 0 1 0 0 0;
              0 0 0 0 1 0];
model(5).D = [0;0;0];

%% Uncertainty of system
n = 6;
n_u = 1;
n_y = 6;
n_mode = 4;
N_A = 0*ones(n,n,n_mode);
N_B = 0*ones(n,n_u,n_mode);


for sc = 1:5  % different scenarios
load('4zone_Radiant_systemmat_300s')
clear C D f C_f D_f f_f N_C N_D Nf_C Nf_D
A(:,:,1) = SYS1.A;
A(:,:,2) = SYS2.A;
A(:,:,3) = SYS3.A;
A(:,:,4) = SYS4.A;
for i = 1:4
B(:,:,i) = [0;0;0;0;0;0];
C(:,:,i) = model(sc).C;
n_y = size(model(sc).C,1);
N_C(:,:,i) = zeros(n_y,n);
D(:,:,i) = model(sc).D;
N_D(:,:,i) = 0*model(sc).D;
end
f(:,1) = SYS1.B; 
f(:,2) = SYS2.B; 
f(:,3) = SYS3.B; 
f(:,4) = SYS4.B; 

n_y = size(model(sc).C,1);
g = zeros(n_y,4);
N_g = zeros(n_y,4);

num_m = size(A,3);

N_f = [zeros(2,4);0.1*f(3:6,:)]; %Uncertainty normalization calculated from problem setup


SYS = StateSpace(A,B,C,D,f,g);
N = StateSpace(N_A,N_B,N_C,N_D,N_f,N_g);

%% Fault Model Parameters
load('4zone_Radiant_Faultmat_300s')
n = 6;
n_u = 1;
n_mode = 2;
Nf_A = 0*ones(n,n,n_mode);
Nf_B = 0*ones(n,n_u,n_mode);
A_f(:,:,1) = SYSF1.A;
A_f(:,:,2) = SYSF2.A;
for i = 1:2
B_f(:,:,i) = [0;0;0;0;0;0];
C_f(:,:,i) = model(sc).C;
Nf_C(:,:,i) = zeros(size(model(sc).C));
D_f(:,:,i) = model(sc).D;
Nf_D(:,:,i) = 0*model(sc).D;
end
f_f(:,1) = SYSF1.B; 
f_f(:,2) = SYSF2.B;  
g_f = zeros(n_y,2);
Nf_g = zeros(n_y,2);
Nf_f = [zeros(2,2);0.1*f_f(3:6,:)]; %Uncertainty normalization calculated from problem setup
SYS_f = StateSpace(A_f,B_f,C_f,D_f,f_f,g_f);
Nf = StateSpace(Nf_A,Nf_B,Nf_C,Nf_D,Nf_f,Nf_g);

%% Bounds

pn_bound = [0 0];
input_low = 0;
input_high = 10;
state_low = [15 15];
state_high = [19 19];
solver = 'cplex';
t1 = tic;

for T = 1:7
    switch SWA
        case 'uswa' 
            mn_bound = [0.05 0.05];
            [Decision,sol] = Tdetect_uswa_milp_modified(SYS,SYS_f,T,N,Nf ...
            ,mn_bound,pn_bound,input_low, input_high, state_low, state_high,solver)
        case 'swa'
            mn_bound = [0 0];
            [Decision,sol] = Tdetect_swa_milp(SYS,SYS_f, ...
            T,mn_bound,pn_bound,input_low, input_high, state_low, state_high,solver)
    end
DD{T,sc} = Decision;
end
toc(t1)
end