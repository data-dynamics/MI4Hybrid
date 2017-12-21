% Example for invalidation_pwa_milp function
% Define a PWA system => simulate the system 
% => invalidate the system through input and generated output 
%clear,close all,clc
clear,close all,
%addpath('../lib/')

% System Parameters  
A(:,:,1) = [1 0.095;-25 -2];
A(:,:,2) = [1 0.1;-22.5 -2];
B(:,:,1) = [0; 1];
B(:,:,2) = [0; 1];
C(:,:,1) = [5 2];
C(:,:,2) = [1.5 2];
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
T = 10;

% Generate input & switching sequence
%input = 5*randn(1,T);
%input = [-12.0603    4.3306   -0.9212   -2.4405   -2.1919   -8.7977   -3.2080   -7.8442   -3.6496    1.1727];
input = [7.8254   -8.4667   -2.2470   -0.4215   -9.9600    4.2062   -2.0733    9.5609   -1.9545    2.0459];
% Generate data using the system sys, input and switchseq
[y,x,p_noise,m_noise,switchseq]=pwa_sim(sys,input);
figure(1),plot(x(1,:),x(2,:)),hold on,plot([2,-2],[-1,3],'r'),axis([-2,2,-2,2])
xlabel('x1'),ylabel('y1'),legend('states plot','classification hyperplane')

% Model invalidation
[Decision] = invalidation_pwa_milp2(sys,input,y,0,0,10,100, 'gurobi')
%[Decision] = invalidation_pwa_milp(sys,input,y,0,0,10,100, 'gurobi')
% 
% % what if I give a wrong output sequence?
Decision = invalidation_pwa_milp2(sys,input,y+0.1,0,0,10,100, 'gurobi')
  %Decision = invalidation_pwa_milp(sys,input,y+0.1,0,0,10,100, 'gurobi')
% 
% %% include noise in the system
% [y,x,p_noise,m_noise,switchseq]=pwa_sim(sys,input,[],0.1,0.1);
% figure(1),plot(x(1,:),x(2,:)),hold on,plot([2,-2],[-1,3],'r'),axis([-2,2,-2,2])
% xlabel('x1'),ylabel('y1'),legend('states plot','classification hyperplane')
% 
% % Model invalidation
% Decision = invalidation_pwa_milp2(sys,input,y,0.1,0.1,10000,10000, 'gurobi')
% Decision = invalidation_pwa_milp(sys,input,y,0.1,0.1,10000,10000, 'gurobi')
% 
% % lets change to another system
% A(:,:,1) = [11 0.095;-25 -10];
% sys = PWAModel(A,B,C,D,P,M,f,g);
% Decision = invalidation_pwa_milp2(sys,input,y,0.1,0.1,10000,10000, 'gurobi')
% Decision = invalidation_pwa_milp(sys,input,y,0.1,0.1,10000,10000, 'gurobi')