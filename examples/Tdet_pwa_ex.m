% Example for invalidation_pwa_milp function
% Define a PWA system => simulate the system 
% => invalidate the system through input and generated output 
clear,close all,clc
%addpath('../lib/')

% System Parameters  
A(:,:,1) = [1 0.095;-25 -2];
A(:,:,2) = [1 0.1;-22.5 -2];
B(:,:,1) = [0; 1];
B(:,:,2) = [0; 1];
% C(:,:,1) = [1 2];
% C(:,:,2) = [1 -2];
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

Af(:,:,1) = [11 0.095;-25 -10];;
Af(:,:,2) = [1 0.1;12.5 -2];
Bf(:,:,1) = [0; 1];
Bf(:,:,2) = [0; 1];
Cf(:,:,1) = [5 2];
Cf(:,:,2) = [1.5 2];
% Cf(:,:,1) = [1 2];
% Cf(:,:,2) = [1 -2];
Df(:,:,1) = 0;
Df(:,:,2) = 0;
ff = [1,1;2,-1];
gf = [0,0];
Pf(:,:,1) = [1 1];
Pf(:,:,2) = [-1 -1];
Mf = [-1,1];
% Define a system (sys) in StateSpace class
sysf = PWAModel(Af,Bf,Cf,Df,Pf,Mf,ff,gf);

% Time Horizon
T = 5;

% Generate input & switching sequence
%input = 10*randn(1,T);
input = [-12.0603    4.3306   -0.9212   -2.4405   -2.1919   -8.7977   -3.2080   -7.8442   -3.6496    1.1727];

% % Generate data using the system sys, input and switchseq
% [y,x,p_noise,m_noise,switchseq]=pwa_sim(sys,input);
% figure(1),plot(x(1,:),x(2,:)),hold on,plot([2,-2],[-1,3],'r'),axis([-2,2,-2,2])
% xlabel('x1'),ylabel('y1'),legend('states plot','classification hyperplane')

% T-detectability test
[Decision] = Tdetect_pwa_milp_M(sys, sysf, T,0,0,10,100, 'gurobi')
%[Decision] = invalidation_pwa_milp(sys,input,y,0,0,10000,10000, 'gurobi')