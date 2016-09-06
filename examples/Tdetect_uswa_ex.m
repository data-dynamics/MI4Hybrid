%% Example for Tdetectability of USWA 
% Change parameters of model and fault to obtain different results.
clear, close all

%% system model
A(:,:,1) = 0.1*[1 2;-1 0];
A(:,:,2) = 0.1*[1 -1; 0 1];

B(:,:,1) = [0.1; -0.5];
B(:,:,2) = [0.5; 0.5];

C(:,:,1) = eye(2);
C(:,:,2) = eye(2);

D(:,:,1) = zeros(2,1);
D(:,:,2) = zeros(2,1);

f = [0.5 1;0 1];
g = [0 0;0 0];

SYS = StateSpace(A,B,C,D,f,g);

%% fault model
A_f(:,:,1) = 0.1*[1 2;-1 0];
A_f(:,:,2) = 0.1*[1 -1.5; 0 1];

B_f(:,:,1) = [0.1; -0.5];
B_f(:,:,2) = [0.5; 0.5];

C_f(:,:,1) = eye(2);
C_f(:,:,2) = eye(2);

D_f(:,:,1) = zeros(2,1);
D_f(:,:,2) = zeros(2,1);

f_f = [0.8 1;0 1.3];
g_f = [0 0;0 0];

SYS_f = StateSpace(A_f,B_f,C_f,D_f,f_f,g_f);

%% Uncertainty of system
n = size(A,1);
n_y = size(C,1);
n_u = size(B,2);
n_mode = size(A,3);
N_A = 0*ones(n,n,n_mode);
N_B = 0*ones(n,n_u,n_mode);
N_C = 0*ones(n_y,n,n_mode);
N_D = 0*ones(n_y,n_u,n_mode);
N_f = 0.15*ones(n,n_mode);
N_g = 0*ones(n_y,n_mode);
N = StateSpace(N_A,N_B,N_C,N_D,N_f,N_g);

%% Uncertainty of fault
nf_mode = size(A_f,3);
N_A = 0*ones(n,n,nf_mode);
N_B = 0*ones(n,n_u,nf_mode);
N_C = 0*ones(n_y,n,nf_mode);
N_D = 0*ones(n_y,n_u,nf_mode);
N_f = 0.1*ones(n,nf_mode);
N_g = 0*ones(n_y,nf_mode);
Nf = StateSpace(N_A,N_B,N_C,N_D,N_f,N_g);
%% Other parameters
T = 2;
mn_bnd = 1;
pn_bnd = 0.05;
input_low = -2;
input_high = 2;
state_low = [-1000 -1000];
state_high = [1000 1000];
solver = 'cplex';

%% Run T-detectability for uswa
t1 = tic;
[Decision,sol] = Tdetect_uswa_milp(SYS,SYS_f,T,N,Nf ...
    ,mn_bnd,pn_bnd,input_low, input_high, state_low, state_high,solver)
toc(t1)






