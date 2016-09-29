function [Decision,sol,cz,fz,zc_C] = invalidation_swa_milp(SYS,input,output, ...
    mn_bnd,pn_bnd,input_low, input_high, state_low, state_high,solver)

% This function implements MILP-based T-detectability approach for the
% switched state-space models subject to measurement and process noise and
% parameter uncertainty is discussed in "Guaranteed Model-Based Fault
% Detection in Cyber-Physical Systems: A Model Invalidation Approach" by
% F. Harirchi and N. Ozay.
%
% Preleminaries for this function are: YALMIP and (cplex or gurobi).
%
% Inputs/Parameters:
%   SYS -- a user-defined system class (see StateSpace.m)
%   N -- user-defined normalization matrices for system model, which is
%   itself a StateSpace class object(see StateSpace.m)
%   input -- an n_u-by-T matrix specifying the input sequence to the system
%            , where n_u is the number of inputs and T is time horizon.
%   output -- an n_y-by-T matrix specifying the output sequence to the
%             system, where n_y is the number of outputs and T is time horizon.
%   mn_bound -- an n_y-by-1 vector specifying the bound for measurement
%               noise for system, where n_y is the dimension of the output.
%   pn_bound -- an n-by-1 vector specifying the bound for process
%               noise for the system, where n is the dimension of the states.
%   input_high -- an n_u-by-1 vector specifying upperbound on each element
%                  of input, where n_u is the number of inputs.
%   input_low -- an n_u-by-1 vector specifying lowerbound on each element
%                  of input, where n_u is the number of inputs.
%   state_high -- an n-by-1 vector specifying upperbound on each element of
%                  state variables for the system, where n is the number of
%                  states.
%   state_low -- an n-by-1 vector specifying lowerbound on each element of
%                  state variables for the system, where n is the number of
%                  states.
%   solver -- MILP solver to be used to solve the optimization problem
%            -- 'cplex': uses CPLEX solver from IBM (needs to be installed)
%                for more information on cplex see: "http://www-03.ibm.com
%                /software/products/en/ibmilogcpleoptistud"
%            -- 'gurobi': uses Gurobi as the solver, for more information
%                see "http://www.gurobi.com"
% Output:
%   Decision -- the output is one of the following messages:
%               - 'Faulty'
%               - 'Normal'
%   sol -- A solution structure that is the output of YALMIP
%
% Syntax:
%   [Decision,sol] = invalidation_uswa_milp(SYS,N,input,output, ...
%    mn_bnd,pn_bnd,input_low, input_high, state_low, state_high,solver)
%
% Author: MI4Hybrid
% Date: September 26th, 2016


%% Check if the system model is valid for this function
if(strcmp(SYS.mark,'ss')~=1&&strcmp(SYS.mark,'swss')~=1)
    error('The system model must be a state-space model.');
end

%% Input, States, Output dimensions and time horizon
n_y = size(SYS.mode(1).C,1); % output dimension and time horizon
n_u = size(SYS.mode(1).B,2); % input dimension
n_mode = size(SYS.mode,2); % number of system modes
n = size(SYS.mode(1).A,1); % state dimension
T = size(input,2);


%% Set up default values for empty paramters
if(isempty(mn_bnd))
    mn_bnd=zeros(n_y,1);
end
if(isempty(pn_bnd))
    pn_bnd=zeros(n,1);
end
if(isempty(input_low))
    input_low=zeros(n_u,1)-inf;
end
if(isempty(input_high))
    input_high=zeros(n_u,1)+inf;
end
if(isempty(state_low))
    state_low=zeros(n,1)-inf*ones(n,1);
end
if(isempty(state_high))
    state_high=zeros(n,1)+inf*ones(n,1);
end
if(isempty(solver))
    solver='cplex';
end


%% Convert scalars to vectors
if(size(mn_bnd,1)==1&&n_y>1)
    if(size(mn_bnd,2)==1)
        mn_bnd=ones(n_y,1)*mn_bnd;
        warning(['Bound for measurement noise is a scalar, converted to'...
            ' a vector with identical entries.']);
    else
        error(['mn_bnd dimensions are not consistent.']);
    end
end

if(size(pn_bnd,1)==1&&n>1)
    if(size(pn_bnd,2)==1)
        pn_bnd=ones(n,1)*pn_bnd;
        warning(['Bound for process noise is a scalar, converted to'...
            ' a vector with identical entries.']);
    else
        error(['pn_bnd dimensions are not consistent.']);
    end
end

if(size(input_low,1)==1&&n_u>1)
    if(size(input_low,2)==1)
        input_low=ones(n_u,1)*input_low;
        warning(['Lowerbound for input is a scalar, converted to'...
            ' a vector with identical entries.']);
    else
        error(['input_low dimensions are not consistent.']);
    end
end

if(size(input_high,1)==1&&n_u>1)
    if(size(input_high,2)==1)
        input_high=ones(n_u,1)*input_high;
        warning(['Upperbound for input is a scalar, converted to'...
            ' a vector with identical entries.']);
    else
        error(['input_high dimensions are not consistent.']);
    end
end


if(size(state_low,1)==1&&n>1)
    if(size(state_low,2)==1)
        state_low=ones(n,1)*state_low;
        warning(['State lower bound is a scalar, converted to a vector with '...
            'identical entries.']);
    else
        error(['state_low dimensions are not consistent.']);
    end
end

if(size(state_high,1)==1&&n>1)
    if(size(state_high,2)==1)
        state_high=ones(n,1)*state_high;
        warning(['State upper bound is a scalar, converted to a vector with '...
            'identical entries.']);
    else
        error(['state_up dimensions are not consistent.']);
    end
end

%% Check input and output dimensions
if (size(input, 1)~=n_u)
    error('The number of inputs is not correct.');
end

if (size(output, 1)~=n_y)
    error('The number of outputs is not correct.');
end

if (size(output,2)~=T)
    error('The length of output and input sequences must be the same.');
end

%% Check the bounds
if(size(mn_bnd,1)~=n_y||size(mn_bnd,2)~=1)
    error('The number of bounds for measurement noise is not correct.');
end
if(size(pn_bnd,1)~=n||size(pn_bnd,2)~=1)
    error('The number of bounds for process noise is not correct.');
end
if(size(input_low,1)~=n_u||~isvector(input_low))
    error('The number of lowerbounds for inputs is not correct.');
end
if(size(input_high,1)~=n_u||~isvector(input_high))
    error('The number of upperbounds for inputs is not correct.');
end
if(size(state_low,1)~=n||size(state_low,2)~=1)
    error('The number of lowerbounds for states is not correct.');
end
if(size(state_high,1)~=n||size(state_high,2)~=1)
    error('The number of upperbounds for states is not correct.');
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

M_u = state_high(:,1);
M_l = state_low(:,1);
eps = mn_bnd(:,1);
PNB = pn_bnd(:,1);
U_u = input_high(:,1);
U_l = input_low(:,1);

%% Check if the input satisfies admissible bounds
count = 0;
for t = 1:T
    for i = 1:n_u
        if ((input(i,t)<=U_l(i))||(input(i,t)>=U_u(i)))
            count = count+1;
        end
    end
end
if (count ~=0)
    error('The given input does not satisfy the admissible bounds')
end

%% Optimization variables
% binary variables
d = binvar(n_mode,T,'full');
% real variables
cz = sdpvar(n,n_mode,T,'full');  % current state
fz = sdpvar(n,n_mode,T,'full');  % future state
theta = sdpvar(n_y,n_mode,T,'full'); % measurement noise (system)
eta = sdpvar(n,n_mode,T,'full'); % process noise (system)

%% Creating Constraints
Constraint = [];

%% state and output equation constraints:
for t = 2:T-1    % time index
    for sys_ind = 1: n_mode  % system mode index
        Constraint = [Constraint, fz(:,sys_ind,t)-Mode(sys_ind).A*cz(:,sys_ind,t)==d(sys_ind,t)*Mode(sys_ind).B*input(:,t)+eta(:,sys_ind,t)+d(sys_ind,t)*Mode(sys_ind).f];
        
        Constraint = [Constraint, Mode(sys_ind).C*fz(:,sys_ind,t)+d(sys_ind,t)*Mode(sys_ind).D*input(:,t)+d(sys_ind,t)*Mode(sys_ind).g+theta(:,sys_ind,t)==d(sys_ind,t)*output(:,t)];
    end
end

%% integer constraints
for t = 1:T
    Constraint = [Constraint sum(d(:,t))==1];
end

%% current and future state constraints
for t = 1:T-1
    Constraint = [Constraint sum(cz(:,:,t+1),2)==sum(fz(:,:,t),2)];
end

%% admissible set constraints
for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for i = 1:n_y
            Constraint = [Constraint  norm(theta(i,sys_ind,t),inf)<=...
                d(sys_ind,t)*eps(i)];
        end
        for i = 1:n
            Constraint = [Constraint  norm(eta(i,sys_ind,t),inf)<=...
                d(sys_ind,t)*PNB(i)];
        end
    end
end

for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for i = 1:n
            Constraint = [Constraint d(sys_ind,t)* (M_l(i))<=cz(i,sys_ind,t)];
            Constraint = [Constraint cz(i,sys_ind,t)<=d(sys_ind,t)*M_u(i)];
            Constraint = [Constraint d(sys_ind,t)* (M_l(i))<=fz(i,sys_ind,t)];
            Constraint = [Constraint fz(i,sys_ind,t)<=d(sys_ind,t)*M_u(i)];
        end
    end
end

%% Setting up Options
options = sdpsettings('verbose',0,'solver',solver);

%% Solve the problem
sol = optimize(Constraint,[],options);

%% Print Results
if strcmp(solver,'cplex')
    if strcmp (sol.info , 'Infeasible problem (CPLEX-IBM)')
        Decision = ['The model is invalidated (Faulty).'];
    else
        Decision = ['Normal'];
    end
end
if strcmp(solver, 'gurobi')
    if strcmp(sol.info, 'Infeasible problem (GUROBI-GUROBI)')
        Decision = ['Faulty'];
    else
        Decision = ['The model is not invalidated (Normal).'];
    end
end