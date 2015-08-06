function [Decision,sol,Constraint,e,s] = invalidation_sarx_milp(SYS,input,output,...
    mn_bound,input_bound,solver)

% This function implements MILP-based model invalidation approach for the
% switched ARX models discussed in "Model Invalidation for Switched Affine
% Systems with Applications to Fault and Anomaly Detection" in ADHS15 by
% F. Harirchi and N. Ozay.
%
% Preleminaries for this function are: YALMIP and (cplex or gurobi).
%
% Inputs/Parameters:
%   SYS -- a user-defined system class (see ARXmodel.m)
%   input -- input sequence (an n_i-by-T matrix) where n_i is the dimension
%            of input and T is time horizon
%   output -- output sequence (an n_y-by-T matrix) where n_y is the
%             dimension of output and T is time horizon
%   input_bound -- an n_i-by-1 vector specifying p-norm bound on input
%                  values where n_i is number of inputs.
%   mn_bound -- an n_y-by-1 vector specifying the bound for measurement
%               noise where n_y is the dimension of output.
%   solver -- MILP solver to be used to solve the optimization problem.
%           -- 'cplex': uses CPLEX solver from IBM (needs to be installed)
%               for more information on cplex see: "http://www-03.ibm.com
%               /software/products/en/ibmilogcpleoptistud"
%           -- 'gurobi': uses Gurobi as the solver, for more information
%               see "http://www.gurobi.com"
% Output:
%   Decision -- the output is one of the following messages:
%               - 'Input bounds are not satisfied'
%               - 'The model is invalidated'
%               - 'The model is not invalidated'
%   sol -- A solution structure that is the output of YALMIP
%
% Syntax:
%   Decision = invalidation_sarx_milp(SYS,input,output,mn_bound,input_bound);
%   Decision = invalidation_sarx_milp(SYS,input,output,mn_bound,input_bound,solver);
%
% Author: MI4Hybrid
% Date: July 15th, 2015

%% Check if the system model is valid for this function
if(strcmp(SYS.mark,'arx')~=1&&strcmp(SYS.mark,'swarx')~=1)
    error('The system model must be an ARX model.');
end

%% Check if the input, output, and the model are consistent
if(length(input)~=length(output))
    error('The input length and output length are not the same.');
end
if(size(input,1)~=size(SYS.mode(1).C,2))
    error('The input is not consistent with the model.');
end
if(size(output,1)~=size(SYS.mode(1).A,1))
    error('The output is not consistent with the model.');
end

%% Input, States, Output dimensions and time horizon 
[n_y,T] = size(output); % output dimension and time horizon
n_i = size(input,1); % input dimension
n_mode = size(SYS.mode,2); % number of modes

%% Use the default solver if it is not specified
if(nargin==5)
    solver='cplex';
end

%% Set up default values for empty paramters
if(isempty(mn_bound))
    mn_bound=zeros(n_y,1);
end
if(isempty(input_bound))
    input_bound=zeros(n_i,1)+inf;
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

%% Check the bounds
if(length(mn_bound)~=n_y||~isvector(mn_bound))
    error('The number of bounds for measurement noise is not correct.');
end
if(length(input_bound)~=n_i||~isvector(input_bound))
    error('The number of bounds for inputs is not correct.');
end

%% Initiate system modes
for i = 1: n_mode
    Mode(i).A = SYS.mode(i).A;
    Mode(i).C = SYS.mode(i).C;
    Mode(i).f = SYS.mode(i).f;
end
eps = mn_bound;
U = input_bound;

%% Checking input bounds
for i = 1:n_i
    u_norm(1,i) = norm(input(i,:),SYS.input_norm(i));
end
if sum(u_norm > U) >0
    Decision = 'Input bounds are not satisfied';
    return
end

%% Building the constraints

% Building vectors for output equations
n_a = size(Mode(1).A,3);
n_c = size(Mode(1).C,3);
n = max(n_a,n_c);
for i = 1:n_mode
    MODE(i).A(:,:,1) = -eye(n_y);
    for k = 1:n_a
        MODE(i).A(:,:,k+1) = Mode(i).A(:,:,k);
    end
end
for i = 1:n_mode
    sys(i).h = [];
    for j = 1: n_y
        for k = 1:n_a+1
            sys(i).h = [sys(i).h MODE(i).A(:,j,k)];
        end
    end
    sys(i).g = [];
    for j = 1:n_i
        for k = 1:n_c
            sys(i).g = [sys(i).g Mode(i).C(:,j,k)];
        end
    end
    sys(i).f = Mode(i).f;
end

% Define variables
s = binvar(T,n_mode);
e = sdpvar(n_y*(n_a+1),n_mode*T);

Constraint = [];
for t = n+1:T-1   % time index
    for sys_ind = 1: n_mode  % mode index
        for j = 1:n_y   % constraints for output equation
            Constraint = [Constraint -sys(sys_ind).h(j,:)*e(:,(t-1)*n_mode...
                +sys_ind)+s(t,sys_ind)*(sys(sys_ind).g(j,:)*reshape(input(:,...
                t-1:-1:t-n_c)',1,[])'+sys(sys_ind).h(j,:)*reshape(output(:,...
                t:-1:t-n_a)',1,[])')+ sys(sys_ind).f(j) ==0];
        end
        for i=1:n_y
            %constraints for noise variables
            Constraint = [Constraint norm(e((i-1)*(n_a+1)+1:i*(n_a+1),...
                (t-1)*n_mode+sys_ind),inf)<= eps(n_y)*s(t,sys_ind)];
        end
    end
    Constraint = [Constraint sum(s(t,:))==1];
    % constraints for noise variables
    for i=1:n_y
        Constraint = [Constraint sum(e((i-1)*(n_a+1)+2:(i)*(n_a+1),...
            (t-1)*n_mode+1:t*n_mode),2)==sum(e((i-1)*(n_a+1)+1:i*(n_a+1)-1,...
            (t-2)*n_mode+1:(t-1)*n_mode),2)];
    end
end              

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