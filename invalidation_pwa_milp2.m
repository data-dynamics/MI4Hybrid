function [Decision,sol] = invalidation_pwa_milp2(SYS,input,output,pn_bound,mn_bound,...
    input_bound,state_bound,solver)

% This function implements MILP-based model invalidation approach for the
% piecewise switch affine(PWA) system.
%
% Preleminaries for this function are: YALMIP and (cplex or gurobi).
%
% Inputs/Parameters:
%   SYS -- a user-defined system class (see StateSpace.m)
%   input -- input sequence (an n_i-by-T matrix) where n_i is the dimension
%            of input and T is time horizon
%   output -- output sequence (an n_y-by-T matrix) where n_y is the
%             dimension of output and T is time horizon
%   pn_bound -- an n-by-1 vector specifying the bound for processing
%               noise where n is the dimension of states
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
%   Decision = invalidation_pwa_milp2(SYS,input,output,pn_bound,mn_bound,input_bound,state_bound);
%   Decision = invalidation_pwa_milp2(SYS,input,output,pn_bound,mn_bound,input_bound,...
%                       state_bound,solver);
%
% Author: MI4Hybrid
% Modified: Supratim Ghosh
% Date: December 21st, 2017

%% Check if the system model is valid for this function
if(strcmp(SYS.mark,'pwa')~=1)
    error('The system model must be a piecewise switch affine(PWA) model.');
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

%% Input, States, Output dimensions ,hyperplane and time horizon
[n_y,T] = size(output); % output dimension and time horizon
n_i = size(input,1); % input dimension
n_mode = size(SYS.mode,2); % number of modes
n = size(SYS.mode(1).A,1); % state dimension
n_hyper = size(SYS.mode(1).P,1);%number of hyperplane
%% Use the default solver if it is not specified
if(nargin==7)
    solver='gurobi';
end

%% Set up default values for empty paramters
if(isempty(pn_bound))
    pn_bound=zeros(n,1);
end
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
    solver='gurobi';
end

%% Convert scalars to vectors
if(length(pn_bound)==1&&n>1)
    pn_bound=ones(n,1)*mn_bound;
    warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
end
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
if(length(pn_bound)~=n||~isvector(pn_bound))
    error('The number of bounds for processing noise is not correct.');
end
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
    Mode(i).P = SYS.mode(i).P;
    Mode(i).M = SYS.mode(i).M;
end
S = state_bound;
eps = mn_bound;
gamma = pn_bound;

%% Checking input bounds
for i = 1:n_i
    u_norm(1,i) = norm(input(i,:),SYS.input_norm(i));
end
if sum(u_norm > input_bound) >0
    Decision = 'Input bounds are not satisfied';
    return
end

%% Building the constraints

a = binvar(n_mode,T,'full'); %sequence matrix

% real variables
x = sdpvar(n,T,'full');  % current state
eta = sdpvar(n_y,T,'full'); % measurement noise (system)
nu = sdpvar(n,T,'full'); % process noise (system)

% additional variables in case alternate formulation is used
% z = sdpvar(n, n_mode, T, 'full');
% o = sdpvar(n_y,n_mode, T,'full');

%% Creating Constraints
Constraint = [];

% bounds for big - M formulation
M1 = 1000*ones(n_hyper,1); % bound for hyper-plane constraints
M2 = 10000*ones(n,1); % bound for state constraints
M3 = 10000*ones(n_y,1); % bound for output constaints

%% mode, state and output equation constraints:
for t = 1:T-1    % time index
    for sys_ind = 1: n_mode  % system mode index

        Constraint = [Constraint -M2*(1-a(sys_ind,t))<=x(:,t+1)-Mode(sys_ind).A*x(:,t)-Mode(sys_ind).B*input(:,t)-nu(:,t)-Mode(sys_ind).f];
        Constraint = [Constraint M2*(1-a(sys_ind,t))>=x(:,t+1)-Mode(sys_ind).A*x(:,t)-Mode(sys_ind).B*input(:,t)-nu(:,t)-Mode(sys_ind).f];
        
        Constraint = [Constraint -M3*(1-a(sys_ind,t))<=Mode(sys_ind).C*x(:,t)+Mode(sys_ind).D*input(:,t)+Mode(sys_ind).g+eta(:,t)-output(:,t)];
        Constraint = [Constraint M3*(1-a(sys_ind,t))>=Mode(sys_ind).C*x(:,t)+Mode(sys_ind).D*input(:,t)+Mode(sys_ind).g+eta(:,t)-output(:,t)];

%     Alternate formulation using Bemporad and Morari's paper

%         Constraint = [Constraint, x(:,t+1) == sum(z(:,:,t),2)];
% %        for b = 1:n_mode
%         Constraint = [Constraint, z(:,sys_ind,t)<=M2*(a(sys_ind,t))];
%         Constraint = [Constraint, z(:,sys_ind,t)>= -M2*(a(sys_ind,t))];
%         Constraint = [Constraint, z(:,sys_ind,t)-Mode(sys_ind).A*x(:,t)-Mode(sys_ind).B*input(:,t)-nu(:,t)-Mode(sys_ind).f<= M2*(1-a(sys_ind,t))];
%         Constraint = [Constraint, z(:,sys_ind,t)-Mode(sys_ind).A*x(:,t)-Mode(sys_ind).B*input(:,t)-nu(:,t)-Mode(sys_ind).f>= -M2*(1-a(sys_ind,t))];
%         Constraint = [Constraint, output(:,t) == sum(o(:,:,t),2)];
%         Constraint = [Constraint, o(:,sys_ind,t)<=M3*(a(sys_ind,t))];
%         Constraint = [Constraint, o(:,sys_ind,t)>= -M3*(a(sys_ind,t))];
%         Constraint = [Constraint, o(:,sys_ind,t)-Mode(sys_ind).C*x(:,t)-Mode(sys_ind).D*input(:,t)-eta(:,t)-Mode(sys_ind).g<= M2*(1-a(sys_ind,t))];
%         Constraint = [Constraint, o(:,sys_ind,t)-Mode(sys_ind).C*x(:,t)-Mode(sys_ind).D*input(:,t)-eta(:,t)-Mode(sys_ind).g>= -M2*(1-a(sys_ind,t))];
%         value(a)
%         disp('iteration')
%         disp(t)
    end
end

%% integer constraints
for t = 1:T
    Constraint = [Constraint sum(a(:,t))==1];
end

% last output constraints

  Constraint = [Constraint -M3*(1-a(sys_ind,T))<=Mode(sys_ind).C*x(:,T)+Mode(sys_ind).D*input(:,T)+Mode(sys_ind).g+eta(:,T)-output(:,T)];
  Constraint = [Constraint M3*(1-a(sys_ind,T))>=Mode(sys_ind).C*x(:,T)+Mode(sys_ind).D*input(:,T)+Mode(sys_ind).g+eta(:,T)-output(:,T)];

%  Alternate formulation using Bemporad and Morari's paper

%    Constraint = [Constraint, output(:,T) == sum(o(:,:,T),2)];
%    Constraint = [Constraint, o(:,sys_ind,T)<=M3*(a(sys_ind,T))];
%    Constraint = [Constraint, o(:,sys_ind,T)>= -M3*(a(sys_ind,T))];
%    Constraint = [Constraint, o(:,sys_ind,T)-Mode(sys_ind).C*x(:,T)-Mode(sys_ind).D*input(:,T)-eta(:,T)-Mode(sys_ind).g<= M2*(1-a(sys_ind,T))];
%    Constraint = [Constraint, o(:,sys_ind,T)-Mode(sys_ind).C*x(:,T)-Mode(sys_ind).D*input(:,T)-eta(:,T)-Mode(sys_ind).g>= -M2*(1-a(sys_ind,T))];

for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        Constraint = [Constraint, Mode(sys_ind).P*x(:,t)+Mode(sys_ind).M<= M1*(1-a(sys_ind,t))];
    end
end
        
%% admissible set constraints
for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for i = 1:n_y
            Constraint = [Constraint  norm(eta(i,t),inf)<=eps(i)];
        end
        for i = 1:n
            Constraint = [Constraint  norm(nu(i,t),inf)<=gamma(i)];
        end
    end
end

for t = 1:T    % time index
    for i = 1:n
        Constraint = [Constraint, (-S(i))<= x(i,t)];
        Constraint = [Constraint, x(i,t)<= S(i)];
    end
end

% Use this only if the alternate formulation is used
% for t = 1:T    % time index
%         for i = 1:n
%             Constraint = [Constraint -S(i)<=z(:,:,t)];
%             Constraint = [Constraint z(:,:,t)<= S(i)];
%         end
% end

%% Setting up Options
options = sdpsettings('verbose',0,'solver',solver);
options.gurobi.Heuristics=0.7;


%% Solve the problem
sol = optimize(Constraint,[],options);
value(sol)
AD = value(a)
XN = value(x)

%% Print Results
if strcmp(solver,'cplex')
    if strcmp (sol.info , 'Infeasible problem (CPLEX-IBM)')
        Decision = ['The model is invalidated.'];
        flag = 0;
    else
        Decision = ['The model is not invalidated'];
        flag = 1;
    end
end
if strcmp(solver, 'gurobi')
    if strcmp(sol.info, 'Infeasible problem (GUROBI-GUROBI)')
        Decision = ['The model is invalidated'];
        flag = 0;
    else
        Decision = ['The model is not invalidated.'];
        flag = 1;
    end
end