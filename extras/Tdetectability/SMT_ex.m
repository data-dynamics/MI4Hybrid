%% SMT example

%% Clear
clear, close all

addpath('../../lib/')

%% System model

A(:,:,1) = [0.8 1;0 0.1];
A(:,:,2) = [0.1 0.3; 1 0];
B(:,:,1) = [1;1];
B(:,:,2) = [1;1];
C(:,:,1) = eye(2);
C(:,:,2) = eye(2);
D(:,:,1) = [0; 0];
D(:,:,2) = [0; 0];
F(:,1) = [3;3];
F(:,2) = [1;1];
nbnd = [-0.1 0.1];
xbnd = [0 inf];

sys = StateSpace(A,B,C,D,F,[]);


%% Fault model 

A_f(:,:,1) = [0.8 1;0 0.1];
B_f(:,:,1) = [1;1];
C_f(:,:,1) = eye(2);
D_f(:,:,1) = [0; 0];
F_f(:,1) = [10;10];


sys_f = StateSpace(A_f,B_f,C_f,D_f,F_f,[]);

nfbnd = [-0.1 0.1];
xfbnd = [0 inf];
T = 3;


filename= 'trial1.smt2';

SMTgen(filename,T,sys,xbnd,nbnd,sys_f,xfbnd,nfbnd)