%% Example for Model Invalidation of USWA 
% Change parameters of model and fault to obtain different results.
clear, close all
addpath('../lib/')

%% system model
A(:,:,1) = 0.3*[1 2 1;-1 0 1;0 1 0];
A(:,:,2) = 0.3*[1 -1 0;0 0 1; 1 0 1];

B(:,:,1) = [0.1; -0.5;0];
B(:,:,2) = [0.5; 0.5;1];

C(:,:,1) = eye(3);
C(:,:,2) = eye(3);

D(:,:,1) = zeros(3,1);
D(:,:,2) = zeros(3,1);

f = [0.5 1;0 1;0 0];
g = [0 0;0 0; 0 0];

SYS = StateSpace(A,B,C,D,f,g);

%% Uncertainty of system
n = size(SYS.mode(1).A,1);
n_y = size(SYS.mode(1).C,1);
n_u = size(SYS.mode(1).B,2);
n_mode = size(SYS.mode,2);
N_A = 0.2*ones(n,n,n_mode);
N_B = 0*ones(n,n_u,n_mode);
N_C = 0*ones(n_y,n,n_mode);
N_D = 0*ones(n_y,n_u,n_mode);
N_f = 0*ones(n,n_mode);
N_g = 0*ones(n_y,n_mode);
N = StateSpace(N_A,N_B,N_C,N_D,N_f,N_g);

eps_p = 0.1;
eps_m = 0.1;
%% Other parameters
mn_bnd = eps_m;
pn_bnd = eps_p;
input_low = -2;
input_high = 2;
state_low = -1000;
state_high = 1000;
solver = 'cplex';


% Generate data from system
x = ones(n,1);
input = (rand(1,10)-0.5)*2*input_high;
modeseq = randi(2,1,10);
Status = 'Healthy';
switch Status
    case 'Faulty'
        eta = (rand(n,10)-0.5)*2*(eps_p+0.15);
    case 'Healthy'
        eta = (rand(n,10)-0.5)*2*(eps_p);
end
theta = (rand(n,10)-0.5)*2*eps_m;
for i = 1:10
    x = (SYS.mode(modeseq(i)).A+N.mode(modeseq(i)).A.*((rand(n,n)-0.5)*2))*x+SYS.mode(modeseq(i)).B*input(i)+SYS.mode(modeseq(i)).f+eta(:,i); 
    output(:,i) = SYS.mode(modeseq(i)).C*x+SYS.mode(modeseq(i)).D*input(i)+SYS.mode(modeseq(i)).g+theta(:,i); 
end
state_bound = [state_low;state_high];
%% Run Model Invalidation for uswa
[Decision,sol] = invalidation_uswa_milp(SYS,N,input,output, ...
    mn_bnd,pn_bnd,input_low, input_high, state_low, state_high,solver);
Decision

[Decision,sol] = invalidation_swa_milp(SYS,input,output, ...
    mn_bnd,pn_bnd,input_low, input_high, state_low, state_high,solver);
Decision






