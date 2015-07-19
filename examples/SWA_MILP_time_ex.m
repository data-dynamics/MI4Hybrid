%% Timing for SWA_MILP function
clear,close all
addpath('../lib/')

% System Parameters
A(:,:,1) = [1 0.095;-25 -2];
A(:,:,2) = [1 0.1;-22.5 -2];
B(:,:,1) = [0; 1];
B(:,:,2) = [0; 1];
C(:,:,1) = [1 0];
C(:,:,2) = [1 0];
D(:,:,1) = 0;
D(:,:,2) = 0;
eps = 0.3; % validates the model
% eps = 0.4; % invalidates the model

% Define a system (sys) in StateSpace class
sys = StateSpace(A,B,C,D,[0 0;0 0],[0 0],[inf inf],inf,zeros(2,2),1,inf,...
    [inf inf]);

% Time Horizon
for T = 100:100:1000;

for i = 1:30
% Generate input & switching sequence
input = 10*randn(1,T);
switchseq = randi(2,1,T);

% Generate data using the system sys, input and switchseq
[y,p_noise,m_noise,switchseq]=swss_sim(sys,input,[],[],eps,[],[10 100],...
    switchseq,0);

% Model invalidation
t1 = tic;
Decision = SWA_MILP(sys,y,input,inf,1000,[10 100],eps, 'cplex')
t(i,T/100) = toc(t1);
end

end

for j = 1:10

t_mean(j)= mean(t(:,j));
t_var(j) = std(t(:,j));
