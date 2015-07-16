% Example for SS invalidation.

clear all
warning off
%addpath 'C:/Users/Administrator/Documents/MATLAB/Summer/MI4Hybrid2/lib'

% A test model.
A = [1 0.095;-25 -2];
B = [0; 1];
C = [1 0];
D = 0;

%{
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
%}

% Creat a system class.
sys = StateSpace(A,B,C,D,[],[],inf,inf,zeros(2,2),1,inf,inf);

% Creat an input sequence.
T=100;
input = 10*randn(1,T);

pn_bound_test=[0.03 0.01]*0;
mn_bound_test=0.3;
state_bound=[10 100];

%{
p=1;
sys2=StateSpace(A*p,B*p,C*p,D*p,g*p,f*p,pn_norm,mn_norm,Ep,Em,input_norm,state_norm);
%}

% Apply invalidation algorithm with the simulation function.
count=0;
for i=1:20
    [output,~,~]=swss_sim(sys,input,[],0,0.5,1000);
    result=InvalidationSS(sys,input,output,pn_bound_test,mn_bound_test,state_bound);
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


