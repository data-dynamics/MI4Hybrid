%% jet engine model invalidation (Assumption: fixed uncertainty for both
%% system and fault models)

clear all
addpath '../lib'

%% Specify model and fault parameters
Ts=0.2;
m_f=0.1;
sigma = (rand-0.5)*2*0.12; % the uncertainty term for both models

BP=[0 0];
BM=[0.05 0.05];
LS = [-5;-5];
US = [5;5];

degmat=[1 0; 2 0; 3 0; 0 1; 0 0];
coeffmat=[1 -1.5*Ts -0.5*Ts -Ts sigma;
    3*Ts 0 0 1-Ts 0];

degmatf=[1 0; 2 0; 3 0; 0 1; 0 0];
coeffmatf=[1 -1.5*Ts -0.5*Ts -Ts m_f+sigma;
    3*Ts 0 0 1-Ts 0];

sys=PolyModel(degmat,coeffmat);
sysf=PolyModel(degmatf,coeffmatf);

rng(666);

%% Generate mixed data (Fault occurs at t=91).
[output,~,~]=poly_sim(sys,randn(0,90),[2;3],BP,BM,[],US,0);
outputf(:,1:90)=output;

new_state=output(:,90);
[output,~,~]=poly_sim(sysf,randn(0,11),new_state,BP,BM,[],US,0);
outputf(:,91:100)=output(:,2:11);

%% Model invalidation in receding horizon manner

% initialization
count1=0;
point_rec=[];  % faulty sample times
yalmip('clear');
count = 1;
T=7; % time horizon size
% model invalidation
for t=7:100
    y=outputf(:,t-T+1:t);
    u = zeros(0,7);
    tic
    [miflag, EPS]=invalidation_poly(sys,u,y,BP, BM, LS, US);
    t
    EPS
    time(count)=toc;
    if(miflag==1)
        fprintf('Invalidated');
        count1=count1+1;
        point_rec=[point_rec,t];
    end
    count = count + 1;
end

count1
avg_time = mean(time)
