%% SMT example

%% Clear
clear, close all

addpath('../../lib/')

%% System model
model(1).C = eye(6);
model(1).D = [0;0;0;0;0;0];
model(2).C = [0 1 0 0 0 0;
              0 0 1 0 0 0;
              0 0 0 1 0 0;
              0 0 0 0 1 0;
              0 0 0 0 0 1];
model(2).D = [0;0;0;0;0];
model(3).C = [0 0 1 0 0 0;
              0 0 0 1 0 0;
              0 0 0 0 1 0;
              0 0 0 0 0 1;
              0 0 0 0 0 0];
model(3).D = [0;0;0;0;0];
model(4).C = [1 0 0 0 0 0;
              0 1 0 0 0 0;
              0 0 1 0 0 0;
              0 0 0 1 0 0;
              0 0 0 0 0 0];
model(4).D = [0;0;0;0;0];
model(5).C = [0 1 0 0 0 0;
              0 0 1 0 0 0;
              0 0 0 0 1 0];
model(5).D = [0;0;0];
model(6).C = [0 1 0 0 0 0;
              0 0 0 0 1 0];
model(6).D = [0;0];
model(7).C = [0 1 0 0 0 0];
model(7).D = [0];


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
load('4zone_Radiant_Faultmat_300s')
A_f(:,:,1) = SYSF1.A;
A_f(:,:,2) = SYSF2.A;
for i = 1:2
B_f(:,:,i) = [0;0;0;0;0;0];
C_f(:,:,i) = model(kk).C;
D_f(:,:,i) = model(kk).D;
end
f_f(:,1) = SYSF1.B; 
f_f(:,2) = SYSF2.B;  
g_f = zeros(mm,2);

delta_f = [zeros(2,2);0.1*f_f(3:6,:)]; 
SYS_f = StateSpace(A_f,B_f,C_f,D_f,f_f,g_f);

mn_bound = [0.05 0.05];
pn_bound = [0 0];
input_low = 0;
input_up = 10;
state_low = [15 15];
state_up = 19;
un_bound = [0.5 0.5];
solver = 'cplex';

% xfbnd = [15 18.5];
t1 = tic;
for T = 8:8
    T
% [Decision] = Tdetect_uuswa_milp(SYS,SYS_f,delta,delta_f,T,mn_bound, ...
%     pn_bound,input_bound,state_bound, un_bound, solver)
[Decision,sol,z_c,z_f] = Tdetect_swa_milp(SYS,SYS_f, ...
    T,mn_bound,pn_bound,input_low, input_up, state_low, state_up, solver)
DD{T,kk} = Decision;
end
toc(t1)
end