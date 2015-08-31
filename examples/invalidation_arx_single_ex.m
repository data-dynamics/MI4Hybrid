%% Example of Non-Switched ARX Model Invalidation

clear all
addpath '../lib'

% Changing this will change system parameters for the invalidation test
% only. e.g. Setting A_factor=0.9, C_factor=0.9, and f_factor=0.9 means
% A*90%, C*90%, and f*90% for all A matrices, C matrices, and vectors f.
A_factor=1;
C_factor=1;
f_factor=1;

% Define system parameter.
A(:,:,1)=[0.4747 0.0628;-0.3424 1.2250];
A(:,:,2)=[0.5230 -0.1144;0.3574 -0.2513];
C(:,:,1)=[1 0 0;0 0 0];
C(:,:,2)=[1 0 0;0 0 0];
f=[0;0];

% Set up noise parameters to generate noisy data.
pn_norm=[inf inf];
mn_norm=[inf inf];
pn_bound=[0.15 0.15]*0;
mn_bound=[0.5 0.5];

% Creat the system model.
sys=ARXmodel(A,C,f,pn_norm,mn_norm);
sys2=ARXmodel(A*A_factor,C*C_factor,f*f_factor,pn_norm,mn_norm);

% Creat random input sequence using normal distribution.
T=100;
degree=size(sys.mode.A,3);
n_in=size(sys.mode.C,2); % Input dimension.
input=[zeros(n_in,degree) randn(n_in,T-degree)];

% Run the example for 20 times.
for i=1:20

    % Run simulation to obtain I/O data.
    [output,~,~,~]=swarx_sim(sys,input,[],pn_bound,mn_bound);

    % Apply the invalidation function.
    % The noise bound here is smaller than the bound used to generate data.
    result=invalidation_arx(sys2,input,output,pn_bound,mn_bound*0.95,'mosek');

    % Display invalidation result.
    if(result==1)
        disp('Validated');
    elseif(result==0)
        disp('Invalidated');
    end

end