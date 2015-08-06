function [Decision,sol] = invalidation_swa_milp(SYS,input,output,mn_bound,...
    input_bound,state_bound,solver)

% This function implements MILP-based model invalidation approach for the
% switched state-space models discussed in "Model Invalidation for Switched
% Affine Systems with Applications to Fault and Anomaly Detection" by
% F. Harirchi and N. Ozay.
%
% Preleminaries for this function are: YALMIP and (cplex or gurobi).
%
% Inputs/Parameters:
%   SYS -- a user-defined system class (see StateSpace.m)
%   input -- input sequence (an n_i-by-T matrix) where n_i is the dimension
%            of input and T is time horizon
%   output -- output sequence (an n_y-by-T matrix) where n_y is the
%             dimension of output and T is time horizon
%   mn_bound -- an n_y-by-1 vector specifying the bound for measurement
%               noise where n_y is the dimension of output
%   input_bound -- an n_i-by-1 vector specifying p-norm bound on input
%                  values where n_i is number of inputs
%   state_bound -- an n-by-1 vector specifying infinity norm bound on state
%                  variables, where n is number of states
%   solver -- MILP solver to be used to solve the optimization problem
%            -- 'cplex': uses CPLEX solver from IBM (needs to be installed)
%                for more information on cplex see: "http://www-03.ibm.com
%                /software/products/en/ibmilogcpleoptistud"
%            -- 'gurobi': uses Gurobi as the solver, for more information
%                see "http://www.gurobi.com"
% Output:
%   Decision -- the output is one of the following messages:
%               - 'Input bounds are not satisfied'
%               - 'The model is invalidated'
%               - 'The model is not invalidated'
%   sol -- A solution structure that is the output of YALMIP
%
% Syntax:
%   Decision = invalidation_swa_milp(SYS,input,output,mn_bound,input_bound,state_bound);
%   Decision = invalidation_swa_milp(SYS,input,output,mn_bound,input_bound,...
%                       state_bound,solver);
%
% Author: MI4Hybrid
% Date: July 13th, 2015

%% Check if the system model is valid for this function
if(strcmp(SYS.mark,'ss')~=1&&strcmp(SYS.mark,'swss')~=1)
    error('The system model must be a state-space model.');
end

%% Check if the input, output, and the model are consistent
if(length(input)~=length(output))
    error('The input length and output length are not the same.');
end
if(size(input,1)~=size(SYS.mode(1).B,2))
    error('The input is not consistent with the model.');
end
if(size(output,1)~=size(SYS.mode(1).C,1))
    error('The output is not consistent with the model.');
end

%% Input, States, Output dimensions and time horizon
[n_y,T] = size(output); % output dimension and time horizon
n_i = size(input,1); % input dimension
n_mode = size(SYS.mode,2); % number of modes
n = size(SYS.mode(1).A,1); % state dimension

%% Use the default solver if it is not specified
if(nargin==6)
    solver='cplex';
end

%% Set up default values for empty paramters
if(isempty(mn_bound))
    mn_bound=zeros(n_y,1);
end
if(isempty(input_bound))
    input_bound=zeros(n_i,1)+inf;
end
if(isempty(state_bound))
    state_bound=zeros(n,1)+inf;
end
if(isempty(solver))
    solver='cplex';
end

%% Convert scalars to vectors
if(length(mn_bound)==1&&n_y>1)
    mn_bound=ones(n_y,1)*mn_bound;
    warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
end
if(length(input_bound)==1&&n_i>1)
    input_bound=ones(n_i,1)*input_bound;
    warning(['Input bound is a scalar, converted to a vector with '...
        'identical entries.']);
end
if(length(state_bound)==1&&n>1)
    state_bound=ones(n,1)*state_bound;
    warning(['State bound is a scalar, converted to a vector with '...
        'identical entries.']);
end

%% Check the bounds
if(length(mn_bound)~=n_y||~isvector(mn_bound))
    error('The number of bounds for measurement noise is not correct.');
end
if(length(input_bound)~=n_i||~isvector(input_bound))
    error('The number of bounds for inputs is not correct.');
end
if(length(state_bound)~=n||~isvector(state_bound))
    error('The number of bounds for states is not correct.');
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
M = state_bound;
eps = mn_bound;

%% Checking input bounds
for i = 1:n_i
    u_norm(1,i) = norm(input(i,:),SYS.input_norm(i));
end
if sum(u_norm > input_bound) >0
    Decision = 'Input bounds are not satisfied';
    return
end

%% Building the constraints
for i = 1:n_mode
    sys(i).g = [-Mode(i).C -ones(n_y,1)];
    sys(i).q = [-Mode(i).D ones(n_y,1)];
    sys(i).h = [-Mode(i).A eye(n)];
    sys(i).l = [-Mode(i).B];
end

s = binvar(T,n_mode);
e = sdpvar(n_mode*T,n_y);
X = sdpvar(n_mode*T,2*n);

Constraint = [];
for t = 1:T-1    % time index
    for sys_ind = 1: n_mode  % mode index
        for j = 1:n_y   % constraints for output equation
            Constraint = [Constraint sys(sys_ind).g(j,:)*[X((t-1)* ...
                n_mode+sys_ind,1:n) e((t-1)*n_mode+sys_ind,j)]'+ ...
                s(t,sys_ind)*(sys(sys_ind).q(j,:)* ...
                [input(:,t) output(j,t)]'- Mode(sys_ind).f(j))==0];
        end
        for state_ind = 1:n  % constraints for state equation
            Constraint = [Constraint sys(sys_ind).h(state_ind,:)* ...
                X((t-1)*n_mode+sys_ind,:)'+ ...
                s(t,sys_ind)*(sys(sys_ind).l(state_ind,:)*input(:,t)'- ...
                Mode(sys_ind).g(state_ind)) ==0];
            Constraint = [Constraint norm(X((t-1)*n_mode+sys_ind, ...
                state_ind),inf)<=M(state_ind)*s(t,sys_ind)];
        end
        Constraint = [Constraint norm(e((t-1)*n_mode+sys_ind,:),inf)<=...
            eps*s(t,sys_ind)];
    end
    Constraint = [Constraint sum(s(t,:))==1];
    Constraint = [Constraint sum(X((t-1)*n_mode+1:t*n_mode,n+1:2*n),1)...
        ==sum(X(t*n_mode+1:(t+1)*n_mode,1:n),1)];
end

% Constraint for the last output equation
for sys_ind = 1: n_mode
    for j = 1:n_y   % constraints for output equation
        Constraint = [Constraint sys(sys_ind).g(j,:)*[X((T-1)*n_mode+ ...
            sys_ind,n+1:2*n) e((T-1)*n_mode+sys_ind,j)]'+ ...
            s(T,sys_ind)*(sys(sys_ind).q(j,:)*[input(:,T) output(j,T)]' +...
            Mode(sys_ind).f(j))==0];
    end
    Constraint = [Constraint norm(e((T-1)*n_mode+sys_ind,:),inf)<= ...
        eps*s(T,sys_ind)];
    for state_ind = 1:n
        Constraint = [Constraint norm(X((T-1)*n_mode+sys_ind,n+state_ind),inf)<=...
            M(state_ind)*s(T,sys_ind)];
    end
end
Constraint = [Constraint sum(s(T,:))==1];

%% Setting up Options
options = sdpsettings('verbose',0,'solver',solver);

%% Solve the problem
sol = optimize(Constraint,[],options);

%% Print Results
if strcmp(solver,'cplex')
    if strcmp (sol.info , 'Infeasible problem (CPLEX-IBM)')
        Decision = 'The model is invalidated';
    else
        Decision = 'The model is not invalidated';
    end
end
if strcmp(solver, 'gurobi')
    if strcmp(sol.info, 'Infeasible problem (GUROBI-GUROBI)')
        Decision = 'The model is invalidated';
    else
        Decision = 'The model is not invalidated';
    end
end