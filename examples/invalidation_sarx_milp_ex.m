%% Example for invalidation_sarx_milp function
clear,close all
addpath('../lib/')

% Define system parameter.
A(:,:,1)=[0.4747 0.0628;-0.3424 1.2250];
A(:,:,2)=[0.5230 -0.1144;0.3574 -0.2513];
C(:,:,1)=[1 0 0;0 0 0];
C(:,:,2)=[1 0 0;0 0 0];
f=[0;0];

% Set up noise parameters to generate noisy data.
pn_norm=[inf inf];
mn_norm=[inf inf];
pn_bound=[1.2 0.8]*0;
mn_bound=[0.5 0.5];

% Creat the system model.
sys=ARXmodel(A,C,f,pn_norm,mn_norm);

% Creat input sequence.
T=100;
degree=size(sys.mode.A,3);
n_u=size(sys.mode.C,2); % Input dimension.
input=[zeros(n_u,degree) randn(n_u,T-degree)];

for l = 1:10
% Generate data using the system sys, input and switchseq
[ym,p_noise,m_noise,switchseq]=swarx_sim(sys,input,[],[],mn_bound,[],[],0);

% invalidationARX Approach 
result=invalidation_arx(sys,input,ym,pn_bound,mn_bound*0.98)

% MILP Approach
Decision= invalidation_sarx_milp(sys,input,ym,[0.5 0.5]*0.98,[100 100 100],'cplex')
end