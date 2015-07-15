%% Example for SS_MILP function
clear,close all
addpath('../lib/')

% System Parameters
A(:,:,1) = [1 0.095;-25 -2];
B(:,:,1) = [0; 1];
C(:,:,1) = [1 0];
D(:,:,1) = 0;

% Define a system (sys) in StateSpace class
sys = StateSpace(A,B,C,D,[0;0],[0],[inf inf],inf,[],1,inf,[inf inf]);

% Time Horizon
T = 100;

% Generate input & switching sequence
input = 10*randn(1,T);
switchseq = randi(1,1,T);

% Generate data using the system sys, input and switchseq
[y,p_noise,m_noise,switchseq]=swss_sim(sys,input,[],[],0.3,[],[10 100],...
    switchseq,0);

% Model invalidation
Decision = SWA_MILP(sys,y,input,inf,1000,[10 100],0.3, 'cplex')