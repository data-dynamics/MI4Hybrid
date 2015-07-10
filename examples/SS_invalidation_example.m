% Example for SS invalidation.

clear all
addpath 'C:/Users/Administrator/Documents/MATLAB/Summer/MI4Hybrid2/lib'

% A test model.
A=[0.54 0.21 0.23;0.44 0.27 0.24;0.43 0.21 0.31];
B=[0 1 0]';
C=[1 0 0;0 1 0;0 0 1];
D=[0 0 0]';
g=[0 0 0]';
f=[0 0 0]';

% Define noise parameters
pn_norm=[inf inf inf];
mn_norm=[inf inf inf];
pn_bound=[1 4 2];
mn_bound=[5 4 1];

Ep=eye(3);
Em=eye(3);
input_norm=inf;
input_bound=inf;
state_norm=[inf inf inf];
state_bound=[1 1 1]*20;

% Creat a system class.
sys=StateSpace(A,B,C,D,g,f,pn_norm,mn_norm,Ep,Em,input_norm,input_bound,state_norm,state_bound);

% Creat an input sequence.
T=500;
input=ones(1,T);
ini_cond=[];

pn_bound_test=[1 4 2];
mn_bound_test=[5 4 1];
p=0.8;
sys2=StateSpace(A*p,B*p,C*p,D*p,g*p,f*p,pn_norm,mn_norm,Ep,Em,input_norm,input_bound,state_norm,state_bound);

% Apply invalidation algorithm with the simulation function.
count=0;
for i=1:20
    [output,~,~,~]=simulates(sys,input,T,ini_cond,pn_bound,mn_bound);
    result=InvalidationSS(sys2,input,output,pn_bound_test,mn_bound_test);
    if(result==false)
        disp('invalidated');
        count=count+1;
    else
        disp('validated');
    end
end

count


