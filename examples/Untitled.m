% This is a file running for system invalidation. The system model is an ARX
%   model with a single mode (i.e. only one submodel).

clear all
addpath 'C:/Users/Administrator/Documents/MATLAB/Summer/MI4Hybrid2/lib'

% Change this to change system parameters for the invalidation test.
% This will not change the generated noisy data, but only influence the invalidation test.
% e.g. 10% change of parameters means A_factor=0.9 and C_factor=0.9.
% No need to care about f because f is a zero vector in this example.
A_factor=1;
C_factor=1;

% Define system parameter.
A(:,:,1)=[0.4747 0.0628;-0.3424 1.2250];
A(:,:,2)=[0.5230 -0.1144;0.3574 -0.2513];
C(:,:,1)=[1 0 0;0 0 0];
C(:,:,2)=[1 0 0;0 0 0];
f=[0;0];

% Set up noise bound.
pn_norm=[inf inf];
mn_norm=[inf inf];
pn_bound=[12 8];    % The process noise bound for data generation only
mn_bound=[10 5];    % The measurement noise bound for data generation only

% Creat the system model.
sys=ARXmodel(A,C,f,pn_norm,mn_norm);

degree=size(sys.mode.A,3);
n_in=size(sys.mode.C,2); % Input dimension.
n_y=size(sys.mode.A,1); % Output dimension.

% Set up desired input data and time horizon T (T may be greater than the length of input).
T=500;
input=[zeros(n_in,degree) ones(n_in,T-degree)];

% Set up initial conditions if needed.
ini_cond=[];






count=0;
for ccc=1:20

% Run simulation to obtain I/O data.
[y,~,~,~]=simulates(sys,input,T,ini_cond,pn_bound,mn_bound);
    
    
result=InvalidationARX(sys,input,y,pn_bound,mn_bound);

if(result==1)
    disp('Valid');
elseif(result==0)
    disp('Invalid');
    count=count+1;
end






end
count