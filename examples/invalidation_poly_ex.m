%% Example of Pendulum Model with 4th Order Taylor Series Expansion
% This example can be used to generate Fig.1 and Fig.2 in paper "Model
% (In)validation and Fault Detection for Systems with Polynomial
% State-Space Models".
clear all
addpath '../lib'

% Sampling time.
Ts=0.3;

% Set up the descritized polynomial model
param_ratio=1;
degmat=[3 0; 1 0; 0 1];
coeffmat=[0     1    Ts;
    Ts/6  -Ts   0.9];
d_coeffmat = [ 0 0.01 0.01*Ts;
    0 0.01*Ts 0.009];
% Uncertain polynomial model
sys=UnPolyModel(degmat,coeffmat*param_ratio,d_coeffmat);

% Certain polynomial model
% sys=PolyModel(degmat,coeffmat*param_ratio);

% Set up noise and state bounds.
noise_ratio=1;  % This parameter can be changed to make data invalid
pn_bound=[0;0];
mn_bound=[0.1;0.1];
state_bound=[2;1.7];

T=10; % number of samples to be used
rng(666); % start with the same initial seed

% To record the results and time.
n_trials = 5;
count=0;
sol_rec=zeros(n_trials,1);
time=zeros(n_trials,1);

% Run n_trials tests.
for i=1:n_trials
    i
    yalmip('clear');
    
    % Specify the input.
    u=randn(0,T);
    
    % Simulate the I/O data.
    [y,~,~]=poly_sim(sys,u,[2;0],pn_bound,mn_bound*noise_ratio,[],[],0);
    
    % Apply polynomial model invalidation function.
    verbose=0;
    param.relaxOrder=2;
    param.scalingSW=0;
    tic
    [miflag, EPS]=invalidation_poly(sys,u,y,pn_bound,mn_bound,-state_bound,state_bound, ...
        param,verbose);
    EPS
    time(i)=toc;
    sol_rec(i)=miflag;
    if(miflag==1)
        fprintf('Invalidated \n');
    else
        fprintf('Not Invalidated \n');
    end
    count=count+1;
end