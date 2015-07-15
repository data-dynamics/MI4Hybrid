% Example for SS invalidation.

%clear all
%addpath 'C:/Users/Administrator/Documents/MATLAB/Summer/MI4Hybrid2/lib'

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
pn_bound=[0.02 0.03 0.01]*0;
mn_bound=[0.5 0.5 0.5];

Ep=eye(3);
Em=eye(3);
input_norm=inf;
input_bound=inf;
state_norm=[inf inf inf];
state_bound=[2.7159 3.9529 2.6302]*2;

% Creat a system class.
sys=StateSpace(A,B,C,D,g,f,pn_norm,mn_norm,Ep,Em,input_norm,state_norm);

% Creat an input sequence.
T=500;
%input=randn(1,T);
ini_cond=[];

pn_bound_test=[0.02 0.03 0.01]*0;
mn_bound_test=[0.3 0.3 0.3];
p=1;
sys2=StateSpace(A*p,B*p,C*p,D*p,g*p,f*p,pn_norm,mn_norm,Ep,Em,input_norm,state_norm);

% Apply invalidation algorithm with the simulation function.
count=0;
for i=1:20
    [output,~,~]=ss_sim(sys,input,ini_cond,pn_bound,mn_bound);
    result=InvalidationSS(sys2,input,output,pn_bound_test,mn_bound_test,state_bound);
    %{
    result=true;
    bound1(i)=norm(output(1,:),inf);
    bound2(i)=norm(output(2,:),inf);
    bound3(i)=norm(output(3,:),inf);
    %}
    if(result==false)
        disp('invalidated');
        count=count+1;
    else
        disp('validated');
    end
end

count


