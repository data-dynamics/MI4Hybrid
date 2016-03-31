%% Example of Moore-Greitzer Jet Engine Model
% This example can be used to generate Fig.4 in paper "Model 
% (In)validation and Fault Detection for Systems with Polynomial 
% State-Space Models".
clear all
addpath '../lib'

% Sampling time.
Ts=0.2;
% random uncertainty value for both system and fault
sigma = (rand-0.5)*2*0.12; 

% Set up the descritized polynomial Jet model 
degmat=[ 1 0; 0 1; 2 0;3 0;0 0];
coeffmat=[1     -Ts -3*Ts/2 -Ts/2 sigma;
         3*Ts  1-Ts  0 0 0];
d_coeffmat = [ 0 0 0 0 0;
               0 0 0 0 0]; 
% Uncertain polynomial model
sys=UnPolyModel(degmat,coeffmat,d_coeffmat);

% Set up the faulty descritized polynomial Jet model 
degmatf=[ 1 0; 0 1; 2 0;3 0;0 0];
coeffmatf=[ 1  -Ts -3*Ts/2 -Ts/2 sigma+0.1;
          3*Ts  1-Ts  0 0 0];
d_coeffmatf = [ 0 0 0 0 0;
               0 0 0 0 0]; 
% Uncertain polynomial model
sysf=UnPolyModel(degmatf,coeffmatf,d_coeffmatf); 

T = 7;
BM = 0.05;
BMf = 0.05;
LS = -5;
LSf = -5;
US = 5;
USf = 5;
LU = [];
LUf = [];
UU = [];
UUf = [];
BP = 0;
BPf = 0;

[tdflag,epsilon]=tdet_poly(sys,sysf,T,BM,BMf,LS,LSf,US,USf,...
    LU,LUf,UU,UUf,BP,BPf)


