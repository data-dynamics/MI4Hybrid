function [Decision,sol] = invalidation_uswa_milp_modified(SYS,Delta,input,output,mn_bound,...
    input_bound,state_bound,uncertainty_bound,solver)

% This function implements MILP-based model invalidation approach for the
% switched state-space models. This is customized for the examples
% in the paper "Guaranteed Model-Based Fault Detection in Cyber-Physical 
% Systems: A Model Invalidation Approach".
%
% Syntax:
%   [Decision,sol] = invalidation_uuswa_milp_modified(SYS,Delta,input,output,mn_bound,...
%    input_bound,state_bound,uncertainty_bound,solver)
%
% Author: MI4Hybrid
% Date: September 10th, 2016

%% Check if the system model is valid for this function
if(strcmp(SYS.mark,'ss')~=1&&strcmp(SYS.mark,'swss')~=1)
    error('The system model must be a state-space model.');
end

%% Input, States, Output dimensions and time horizon
n_y = size(SYS.mode(1).C,1); % output dimension and time horizon
n_u = size(SYS.mode(1).B,2); % input dimension
n_mode = size(SYS.mode,2); % number of modes
n = size(SYS.mode(1).A,1); % state dimension

%% Use the default solver if it is not specified
if(nargin==6)
    solver='cplex';
end

if(isempty(solver))
    solver='cplex';
end

%% Initiate system modes
for i = 1: n_mode
    Mode(i).A = SYS.mode(i).A;
    Mode(i).B = SYS.mode(i).B;
    Mode(i).C = SYS.mode(i).C;
    Mode(i).D = SYS.mode(i).D;
    Mode(i).f = SYS.mode(i).f;
    Mode(i).g = SYS.mode(i).g;
end
T = size(input,2);
M_l = state_bound(1,1)*ones(n,1);
M_u = state_bound(2,1)*ones(n,1);
eps = mn_bound*ones(n_y,1);
U_u = input_bound(1,1);
U_l = -input_bound(1,1);
d_l = -uncertainty_bound;
d_u = uncertainty_bound;

%% Variables
d = binvar(n_mode,T,'full');
z_c = sdpvar(n,n_mode,T,'full');  
z_f = sdpvar(n,n_mode,T,'full'); 
theta = sdpvar(n_y,n_mode,T,'full');
delta = sdpvar(n_mode,T,'full');

%% Setup Constraints
Constraint = [];

for t = 1:T-1    % time index
    for sys_ind = 1: n_mode  % system mode index       
        Constraint = [Constraint, z_f(:,sys_ind,t)-Mode(sys_ind).A*z_c(:,sys_ind,t)==...
            d(sys_ind,t)*(Mode(sys_ind).f)+Delta(:,sys_ind)*delta(sys_ind,t)]; 
        Constraint = [Constraint, d(sys_ind,t)*output(:,t)-Mode(sys_ind).C*z_f(:,sys_ind,t)==...
            theta(:,sys_ind,t)];                 
    end
end

% integer constraints
for t = 1:T-1
    Constraint = [Constraint sum(d(:,t))==1 ];   
end

% current and future state constraints
for t = 1:T-1
        Constraint = [Constraint sum(z_c(:,:,t+1),2)==sum(z_f(:,:,t),2)];
end

% admissible set constraints
for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for i = 1: n
            Constraint = [Constraint  d(sys_ind,t)* (M_l(i)) <= z_c(i,sys_ind,t)];
            Constraint = [Constraint z_c(i,sys_ind,t) <= d(sys_ind,t)* M_u(i)];
            Constraint = [Constraint  d(sys_ind,t)* (M_l(i)) <= z_f(i,sys_ind,t)];
            Constraint = [Constraint z_f(i,sys_ind,t) <= d(sys_ind,t)* M_u(i)];
        end
            for i = 1:n_y
                Constraint = [Constraint  d(sys_ind,t)* (-eps(i))<=theta(i,sys_ind,t)];
                Constraint = [Constraint theta(i,sys_ind,t) <= d(sys_ind,t)* eps(i)];
            end
            Constraint = [Constraint  d(sys_ind,t)* (d_l)<=delta(sys_ind,t)];
            Constraint = [Constraint delta(sys_ind,t) <= d(sys_ind,t)* d_u];
    end
end

%% Setting up Options
options = sdpsettings('verbose',0,'solver',solver);

%% Solve the problem
sol = optimize(Constraint,[],options);

%% Print Results
if strcmp(solver,'cplex')
    if strcmp (sol.info , 'Infeasible problem (CPLEX-IBM)')
        Decision = ['The model is invalidated'];
    else
        Decision = ['The model is not invalidated'];
    end
end
if strcmp(solver, 'gurobi')
    if strcmp(sol.info, 'Infeasible problem (GUROBI-GUROBI)')
        Decision = ['The model is invalidated'];
    else
        Decision = ['The model is not invalidated'];
    end
end