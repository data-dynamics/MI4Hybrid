%% Example of Non-Switched State-Space Model Invalidation

clear all
addpath '../lib'

% Set up system parameters.
A=[1 0.095;-25 -2];
B=[0;1];
C=[1 0];
D=0;

% Creat a system class.
sys = StateSpace(A,B,C,D,[],[],[inf inf],inf,[],[],inf,[inf inf]);

% Creat an input sequence.
T=100;
input=10*randn(1,T);

% Set up bounds for invalidation function.
pn_bound_test=[0 0];
mn_bound_test=0;
state_bound=[10 100];

%{
p=1;
sys2=StateSpace(A*p,B*p,C*p,D*p,g*p,f*p,pn_norm,mn_norm,Ep,Em,input_norm,state_norm);
%}

% Apply invalidation algorithm with the simulation function.
for i=1:10
    [output,~,~,~]=swss_sim(sys,input,[],[],0.3);
    result=InvalidationSS(sys,input,output,pn_bound_test,mn_bound_test,state_bound);
    %norm(output,inf)
    if(result==false)
        disp('invalidated');
    else
        disp('validated');
    end
end