%% Model generation function for a two-link robotic arm
% The nonlinear model is first converted into a 4th order linear
% differential equation by linearizing it around 4 different operating
% points giving rise to a piece-wise affine system. The continuous system
% thus arising is then discretized using Euler's discertization.


function [A,B,C,D,P,M,f,g] = model_generation(Iz1, Iz2, m1, m2, l1, l2, r1, r2)
syms theta1 theta2 theta1d theta2d;
syms t1 t2;

% system parameters:
alpha = Iz1 + Iz2 + m1*r1^2 + m2*(l1^2+r2^2);
beta = m2*l1*r2;
delta = Iz2 + m2*r2^2;

% system matrices:
a = [alpha+2*beta*cos(theta2) delta+beta*cos(theta2); delta+beta*cos(theta2) delta];
b = [-beta*sin(theta2)*theta2d -beta*sin(theta2)*(theta1d+theta2d); beta*sin(theta2)*theta1d 0];

% expressions for theta1dd and theta2dd
f=-inv(a)*b*[theta1d;theta2d]+inv(a)*[t1;t2];

% State space matrices for linearized model
Ac = [0 1 0 0; diff(f(1),theta1) diff(f(1),theta1d) diff(f(1),theta2) diff(f(1),theta2d); 0 0 0 1; diff(f(2),theta1) diff(f(2),theta1d) diff(f(2),theta2) diff(f(2),theta2d)];
Bc = [0 0; diff(f(1),t1) diff(f(1),t2); 0 0; diff(f(2),t1) diff(f(2),t2)];
Cc = [1 0 0 0; 0 0 1 0];

% Linearized values of A, B, C

% Around point (pi/2, 1, 0, 0)
theta1 = pi/2; theta1d = 1; theta2 = 0; theta2d = 0; t1 = 0; t2 = 0;
A1 = subs(Ac);
B1 = subs(Bc);
C1 = Cc;

% Around point (0, 0, pi/2, 1)
theta1 = 0; theta1d = 0; theta2 = pi/2; theta2d = 1; t1 = 0; t2 = 0;
A2 = subs(Ac);
B2 = subs(Bc);
C2 = Cc;

% Around point (-pi/2, -1, 0, 0)
theta1 = -pi/2; theta1d = -1; theta2 = 0; theta2d = 0; t1 = 0; t2 = 0;
A3 = subs(Ac);
B3 = subs(Bc);
C3 = Cc;

% Around point (0, 0, -pi/2, -1)
theta1 = 0; theta1d = 0; theta2 = -pi/2; theta2d = -1; t1 = 0; t2 = 0;
A4 = subs(Ac);
B4 = subs(Bc);
C4 = Cc;

% Euler Discretization (Forward)
% Sampling time T
Ts = 0.2;

% Mode 1
% Ad1 = eye(size(A1))+Ts*A1;
% Bd1 = B1*Ts;
Ad1 = expm(double(A1)*Ts);
Bd1 = integral(@(t) expm(double(A1)*t)*double(B1),0,Ts,'ArrayValued',true);
Cd1 = C1;

% Mode 2
% Ad2 = eye(size(A2))+Ts*A2;
% Bd2 = B2*Ts;
Ad2 = expm(double(A2)*Ts);
Bd2 = integral(@(t) expm(double(A2)*t)*double(B2),0,Ts,'ArrayValued',true);
Cd2 = C2;

% Mode 3
% Ad3 = eye(size(A3))+Ts*A3;
% Bd3 = B3*Ts;
Ad3 = expm(double(A3)*Ts);
Bd3 = integral(@(t) expm(double(A3)*t)*double(B3),0,Ts,'ArrayValued',true);
Cd3 = C3;

% Mode 4
% Ad4 = eye(size(A4))+Ts*A4;
% Bd4 = B4*Ts;
Ad4 = expm(double(A4)*Ts);
Bd4 = integral(@(t) expm(double(A4)*t)*double(B4),0,Ts,'ArrayValued',true);
Cd4 = C4;

% System Parameters  
A(:,:,1) = double(Ad1);
A(:,:,2) = double(Ad2);
A(:,:,3) = double(Ad3);
A(:,:,4) = double(Ad4);
B(:,:,1) = double(Bd1);
B(:,:,2) = double(Bd2);
B(:,:,3) = double(Bd3);
B(:,:,4) = double(Bd4);
% C(:,:,1) = [1 2];
% C(:,:,2) = [1 -2];
C(:,:,1) = C1;
C(:,:,2) = C2;
C(:,:,3) = C3;
C(:,:,4) = C4;
D(:,:,1) = zeros(2,2);
D(:,:,2) = zeros(2,2);
D(:,:,3) = zeros(2,2);
D(:,:,4) = zeros(2,2);
%f = [1,1;2,-1];
f = zeros(4,4);
%f = [1 0 -1 0; 1 0 -1 0; 0 1 0 -1; 0 1 0 -1];
g = zeros(2,4);
P(:,:,1) = [1 0 -1 0; -1 0 -1 0];
P(:,:,2) = [1 0 -1 0; 1 0 1 0];
P(:,:,3) = [-1 0 1 0; 1 0 1 0];
P(:,:,4) = [-1 0 1 0; -1 0 -1 0];
M = zeros(2,4);

end

