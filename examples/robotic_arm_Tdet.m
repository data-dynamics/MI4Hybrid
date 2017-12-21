% Fault detection for a robotic arm (two-link planar robot)
% Equation 4.11 on page 165 of Richard Murray's book on
% A Mathematical Introduction to Robotic Manipulation
% Use this file in conjunction with the model generation file to get the PWA model
%
clear,close all,clc

%% Define the system parameters

% Moment of intertia tension
Iz1 = 10;
Iz2 = 10;

% Masses of links
m1 = 10;
m2 = 10;

% Lengths of robotic links
l1 = 10;
l2 = 10;

% distances of the center of mass from link bases
r1 = 5;
r2 = 5;

% Get the nominal system parameters
[A,B,C,D,P,M,f,g] = model_generation(Iz1, Iz2, m1, m2, l1, l2, r1, r2);

% Define a system (sys) in StateSpace class
sys = PWAModel(A,B,C,D,P,M,f,g);

% Time Horizon
T = 6;

%% Uncomment this section if you want to run model invalidation
% % Generate input & switching sequence
% input = .0005*randn(2,T);
% 
% % Generate data using the system sys, input and switchseq
% [y,x,p_noise,m_noise,switchseq]=pwa_sim(sys,input);
% 
% 
% % % Model invalidation
% [Decision] = invalidation_pwa_milp2(sys,input,y,0,0,10,pi, 'gurobi')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Faulty model description
% Moment of intertia tension
Izf1 = 10;
Izf2 = 10;

% Masses of links
mf1 = 10;
mf2 = 2.5;

% Lengths of robotic links
lf1 = 10;
lf2 = 2.5;

% distances of the center of mass from link bases
rf1 = 5;
rf2 = 1.25;

% Get the nominal system parameters
[Af,Bf,Cf,Df,Pf,Mf,ff,gf] = model_generation(Izf1, Izf2, mf1, mf2, lf1, lf2, rf1, rf2);

% Define a system (sys) in StateSpace class
sysf = PWAModel(Af,Bf,Cf,Df,Pf,Mf,ff,gf);


[Decision] = Tdetect_pwa_milp_M(sys, sysf, T,0,0,pi,pi, 'gurobi')
