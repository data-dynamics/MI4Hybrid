%% SMT example

%% Clear
clear, close all

addpath('../../lib/')

%% System model
model(1).C = eye(6);
model(1).D = [0;0;0;0;0;0];


for kk = 1:1
load('4zone_Radiant_systemmat_300s')
clear C D f C_f D_f f_f
A(:,:,1) = SYS1.A;
A(:,:,2) = SYS2.A;
A(:,:,3) = SYS3.A;
A(:,:,4) = SYS4.A;
for i = 1:4
B(:,:,i) = [0;0;0;0;0;0];
C(:,:,i) = model(kk).C;
D(:,:,i) = model(kk).D;
end
f(:,1) = SYS1.B; 
f(:,2) = SYS2.B; 
f(:,3) = SYS3.B; 
f(:,4) = SYS4.B; 

mm = size(model(kk).C,1);
g = zeros(mm,4);

num_m = size(A,3);

delta = [zeros(2,4);0.1*f(3:6,:)]; 

SYS = StateSpace(A,B,C,D,f,g);

%% Fault System Parameters
load('4zone_Radiant_weakFaultmat_300s')
A_f(:,:,1) = SYSWF1.A;
A_f(:,:,2) = SYSWF2.A;
A_f(:,:,3) = SYSWF3.A;
A_f(:,:,4) = SYSWF4.A;
for i = 1:4
B_f(:,:,i) = [0;0;0;0;0;0];
C_f(:,:,i) = model(kk).C;
D_f(:,:,i) = model(kk).D;
end
f_f(:,1) = SYSWF1.B; 
f_f(:,2) = SYSWF2.B;  
f_f(:,3) = SYSWF3.B; 
f_f(:,4) = SYSWF4.B;  
g_f = zeros(mm,4);

delta_f = [zeros(2,4);0.1*f_f(3:6,:)]; 
SYS_f = StateSpace(A_f,B_f,C_f,D_f,f_f,g_f);

mn_bound = [0.05 0.05];
pn_bound = [0 0];
input_bound = [10 10; 0 0];
state_bound = [19 19;15 15];
un_bound = [0.5 0.5];
solver = 'cplex';

% xfbnd = [15 18.5];
t1 = tic;
W = 3;
for T = W+1:10
    T
[Decision] = wTdetect_uuswa_milp(SYS,SYS_f,delta,delta_f,T,W,mn_bound, ...
    pn_bound,input_bound,state_bound, un_bound, solver)
DD{T,kk} = Decision;
end
toc(t1)
end