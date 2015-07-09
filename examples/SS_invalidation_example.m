% Example for SS invalidation.

clear all
addpath 'C:/Users/Administrator/Documents/MATLAB/Summer/MI4Hybrid2/lib'

% A segwway model.
A=[-1 0 0 0;0 -2 0 0;0 0 -3 0;0 0 0 -0.5];
B=[0 0 0 0]';
C=[0 -1 0 0;0 0 1 0];
D=[0 0]';
g=[0 0 0 0]';
f=[0 0]';

% Define noise parameters
pn_norm=[inf inf inf inf];
mn_norm=[inf inf];
pn_bound=[0 0 0 0];
mn_bound=[0.3 0.5];

% Creat a system class.
sys=StateSpace(A,B,C,D,g,f,pn_norm,mn_norm);

% Creat an input sequence.
T=100;
input=ones(1,T)*0.5;
ini_cond=[];

pn_bound_test=[0.2 0.3 0.1 0.4];
mn_bound_test=[0.3 0.5];

% Apply invalidation algorithm with the simulation function.
count=0;
for i=1:5
    [output,~,~,~]=simulates(sys,input,T,ini_cond,pn_bound,mn_bound);
    result=InvalidationSS(sys,input,output,pn_bound_test,mn_bound_test);
    if(result==false)
        disp('invalidated');
        count=count+1;
    end
end

count


