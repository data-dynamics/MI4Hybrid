%% Example for invalidation_swa_milp function
clear,close all
addpath('../lib/')

data = 'valid'

switch data
    case 'valid'
        eps_f = 0.5;
    case 'invalid'
        eps_f = 0.3;
end


% System Parameters
A(:,:,1) = [1 0.095;-25 -2];
A(:,:,2) = [1 0.1;-22.5 -2];
B(:,:,1) = [0; 1];
B(:,:,2) = [0; 1];
C(:,:,1) = [1 0];
C(:,:,2) = [1 0];
D(:,:,1) = 0;
D(:,:,2) = 0;
eps = 0.3; % noise bound


% Define a system (sys) in StateSpace class
sys = StateSpace(A,B,C,D,[0 0;0 0],[0 0],[inf inf],inf,zeros(2,2),1,inf,...
    [inf inf]);

% Time Horizon
T = 100;

% Generate input & switching sequence
input = 10*randn(1,T);
switchseq = randi(2,1,T);

% Generate data using the system sys, input and switchseq
[y,p_noise,m_noise,switchseq]=swss_sim(sys,input,[],[],eps_f,[],[10 100],...
    switchseq,0);

% Model invalidation
Decision = invalidation_swa_milp(sys,y,input,inf,1000,[10 100],eps, 'cplex')


