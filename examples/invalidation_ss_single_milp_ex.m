%% Example for invalidation_swa_milp function
clear,close all
addpath('../lib/')

% System Parameters
A(:,:,1) = [1 0.095;-25 -2];
B(:,:,1) = [0; 1];
C(:,:,1) = [1 0];
D(:,:,1) = 0;

pn_bound_test = [];
mn_bound_test = 0.3;
state_bound = [10 100];
mn_bound = 0.3*1.01;
% Define a system (sys) in StateSpace class
sys = StateSpace(A,B,C,D,[0;0],[0],[inf inf],inf,[],1,inf,[inf inf]);

% Time Horizon
T = 100;

for i = 1:10
% Generate input & switching sequence
input = 10*randn(1,T);
switchseq = randi(1,1,T);

% Generate data using the system sys, input and switchseq
[output,p_noise,m_noise,switchseq]=swss_sim(sys,input,[],pn_bound_test,mn_bound,[],state_bound,...
    switchseq,0);

result=invalidation_ss(sys,input,output,pn_bound_test,mn_bound_test,state_bound)

% Model invalidation
Decision = invalidation_swa_milp(sys,input,output,mn_bound_test,1000,state_bound, 'cplex')
end