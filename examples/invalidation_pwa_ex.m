% Example for invalidation_pwa_milp function
% Define a PWA system => simulate the system 
% => invalidate the system through input and generated output 
clear,close all,clc
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
f = [1,1;2,-1];
g = [0,0];
P(:,:,1) = [1 1];
P(:,:,2) = [-1 -1];
M = [-1,1];
% Define a system (sys) in StateSpace class
sys = PWAModel(A,B,C,D,P,M,f,g);

% Time Horizon
T = 100;

% Generate input & switching sequence
input = 10*randn(1,T);

% Generate data using the system sys, input and switchseq
[y,x,p_noise,m_noise,switchseq]=pwa_sim(sys,input);
figure(1),plot(x(1,:),x(2,:)),hold on,plot([2,-2],[-1,3],'r'),axis([-2,2,-2,2])
xlabel('x1'),ylabel('y1'),legend('states plot','classification hyperplane')

% Model invalidation
Decision = invalidation_pwa_milp(sys,input,y,0,10000,10000, 'gurobi')
% what if I give a wrong output sequence?
Decision = invalidation_pwa_milp(sys,input,y+10,0,10000,10000, 'gurobi')
