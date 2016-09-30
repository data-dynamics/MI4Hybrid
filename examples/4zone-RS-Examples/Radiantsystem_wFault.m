% Radiant System with 2 Pumps 4 Zones

%% CLEAR MEMORY
clear,close all

%% Parameter Initialization

T_a = 10;
K_1 = 1/2.1;
K_2 = 1/2.1;
K_3 = 1/2.2;
K_4 = 1/2.2;
K_r1 = 1/0.125;
K_r2 = 1/0.130;
K_r3 = 1/0.125;
K_r4 = 1/0.130;
K_12 = 1/0.16;
K_13 = 1/0.16;
K_24 = 1/0.16;
K_34 = 1/0.16;
K_w1 = 1/0.05;
K_w2 = 0.5/0.07;
C_1 = 1900;
C_2 = 2100;
C_3 = 2000;
C_4 = 1800;
C_r1 = 3000;
C_r2 = 4000;
T_w1 = 18;
T_w2 = 18;

%% State-Space Model
A0 = [-(K_r1+K_r3)/C_r1 0 K_r1/C_r1 0 K_r3/C_r1 0;
      0 -(K_r2+K_r4)/C_r2 0 K_r2/C_r2 0 K_r4/C_r2;
      K_r1/C_1 0 -(K_r1+K_1+K_12+K_13)/C_1 K_12/C_1 K_13/C_1 0;
      0 K_r2/C_2 K_12/C_2 -(K_r2+K_2+K_12+K_24)/C_2 0 K_24/C_2;
      K_r3/C_3 0 K_13/C_3 0 -(K_r1+K_3+K_13+K_34)/C_3 K_34/C_3;
      0 K_r2/C_4 0 K_24/C_4 K_34/C_4 -(K_r2+K_4+K_24+K_34)/C_4];

M(1).A = A0;
M(3).A = A0;
M(3).A(2,2) = M(3).A(2,2) - K_w2/C_r2;
M(2).A = A0;
M(2).A(1,1) = M(2).A(1,1) -K_w1/C_r1;
M(4).A = A0;
M(4).A(2,2) = M(4).A(2,2) - K_w2/C_r2;
M(4).A(1,1) = M(4).A(1,1) -K_w1/C_r1;

M(1).B = [0 0 K_1*T_a/C_1 K_2*T_a/C_2 K_3*T_a/C_3 K_4*T_a/C_4]';
M(3).B = M(1).B;
M(3).B (2,1) = M(3).B(2,1) + K_w2*T_w2/C_r2;
M(2).B = M(1).B;
M(2).B (1,1) = M(2).B(1,1) + K_w1*T_w1/C_r1;
M(4).B = M(1).B;
M(4).B (2,1) = M(4).B(2,1) + K_w2*T_w2/C_r2;
M(4).B (1,1) = M(4).B(1,1) + K_w1*T_w1/C_r1;

M(1).C = eye(6);
M(2).C = eye(6);
M(3).C = eye(6);
M(4).C = eye(6);

% define correct A matrix for each mode
T_s = 300;

sys1= ss(M(1).A,M(1).B,M(1).C,zeros(6,1));
sys1d = c2d(sys1,T_s);
SYSWF1 = mergesys(sys1d);

sys2= ss(M(2).A,M(2).B,M(2).C,zeros(6,1));
sys2d = c2d(sys2,T_s);
SYSWF2 = mergesys(sys2d);

sys3= ss(M(3).A,M(3).B,M(3).C,zeros(6,1));
sys3d = c2d(sys3,T_s);
SYSWF3 = mergesys(sys3d);

sys4= ss(M(4).A,M(4).B,M(4).C,zeros(6,1));
sys4d = c2d(sys4,T_s);
SYSWF4 = mergesys(sys4d);

% Save weak fault model
save('4zone_Radiant_weakFaultmat_300s.mat','SYSWF1','SYSWF2','SYSWF3','SYSWF4')
     
     
     
     
