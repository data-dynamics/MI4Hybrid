%% Example of Non-Switched State-Space Model Invalidation

clear all
addpath '../lib'

% Changing this will change system parameters for the invalidation test
% only. e.g. Setting A_fac=0.9 means A*90% when creating sys2.
A_fac=1;
B_fac=1;
C_fac=1;
D_fac=1;
g_fac=1;
f_fac=1;

% Set up system parameters.
A=[1 0.095;-25 -2];
B=[0;1];
C=[1 0];
D=0;
g=[0;0];
f=0;

% Creat a system class.
sys=StateSpace(A,B,C,D,g,f,[inf inf],inf,[],[],inf,[inf inf]);
sys2=StateSpace(A*A_fac,B*B_fac,C*C_fac,D*D_fac,g*g_fac,f*f_fac,...
    [inf inf],inf,[],[],inf,[inf inf]);

% Creat an input sequence.
T=100;
input=10*randn(1,T);

% Set up bounds for invalidation function.
input_bound=1000;
pn_bound=[0 0];
mn_bound=0.3;
state_bound=[10 100];

% Run the example for 20 times.
for i=1:20
    
    % Run simulation to obtain I/O data.
    [output,~,~,~]=swss_sim(sys,input,[],pn_bound,mn_bound,input_bound,...
        state_bound);
    
    % Apply the invalidation function.
    % The noise bound here is smaller than the bound used to generate data.
    result=invalidation_ss(sys2,input,output,pn_bound,mn_bound*0.95,...
        state_bound,'cplex');
    
    % Display invalidation result.
    if(result==false)
        disp('Invalidated');
    else
        disp('Validated');
    end
    
end