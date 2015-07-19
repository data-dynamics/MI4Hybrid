%% Example Radiant System Fault Detection (Section 4.2 in paper)

clear,close all
addpath('../lib/')

%% System Parameters

A(:,:,1) = [0.54 0.21 0.23; 0.44 0.27 0.24; 0.43 0.21 0.31];
A(:,:,2) = [0.16 0.11 0.12; 0.22 0.23 0.19; 0.22 0.17 0.27];
B(:,:,1) = [0;0;0];
B(:,:,2) = [0;0;0];
C(:,:,1) = [1 0 0];
C(:,:,2) = [1 0 0];
D(:,:,1) = 0;
D(:,:,2) = 0;
g = zeros(3,3);
F = [0 19.8];

num_m = size(A,3);

%% Fault System Parameters
A_f(:,:,1) = [0.54 0.21 0.23; 0.44 0.27 0.24; 0.43 0.21 0.31];
A_f(:,:,2) = [0.16 0.11 0.12; 0.22 0.23 0.19; 0.22 0.17 0.27];
B_f(:,:,1) = [0;0;0];
B_f(:,:,2) = [0;0;0];
C_f(:,:,1) = [1 0 0];
C_f(:,:,2) = [1 0 0];
D_f(:,:,1) = 0;
D_f(:,:,2) = 0;
g_f = zeros(3,3);
F_f = [0 19.8];

% Define a system (sys) in StateSpace class
sys = StateSpace(A,B,C,D,[],F);

% Define faulty system (sys_f) in StateSpace class
sys_f = StateSpace(A_f,B_f,C_f,D_f,[],F_f);

% Time Horizon
T = 200;

% Generate input & switching sequence
input = 10*randn(1,T);
switchseq = randi(num_m,1,T);
eps = 0.3;

% Generate data using the system sys_f, input and switchseq
[y,p_noise,m_noise,switchseq]=swss_sim(sys_f,input,[],[],eps,[],[100 100 100],...
    switchseq,0);

% Model invalidation
Decision = SWA_MILP(sys,y,input,inf,1000,[10 100],0.3, 'cplex')


