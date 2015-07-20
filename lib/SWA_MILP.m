function [Decision, sol] = SWA_MILP(SYS,y,u,p,u_bound,x_bound,mn_bound, solver)

% This function implements MILP-based model invalidation approach for 
% switched state-space models discussed in "Model Invalidation for Switched 
% Affine Systems with Applications to Fault and Anomaly Detection" by 
% F. Harirchi and N. Ozay.
%
% Preleminaries for this function are: YALMIP and (cplex or gurobi).
%
% Inputs/Parameters:
%   SYS -- a user-defined system class (see StateSpace.m)
%   y -- Output sequence
%   u -- input sequence (an n_u-by-T matrix) where n_u is the dimension
%            of input and T is time horizon
%   p -- a scalar denoting the type of the norm to be used to put 
%         constraints on input.
%   u_bound -- an n_u-by-1 vector specifying p-norm bound on input values, 
%               where n_u is number of inputs.
%   x_bound -- an n-by-1 vector specifying infinity norm bound on state 
%               variables, where n is number of states.
%   mn_bound -- an n_y-by-1 vector specifying the bound for process noise
%                where n_y is the dimension of output.
%   solver -- MILP solver to be used to solve the optimization problem.
%            -- 'cplex': uses CPLEX solver from IBM (needs to be installed)
%                 for more information on cplex see: "http://www-03.ibm.com
%                 /software/products/en/ibmilogcpleoptistud"
%            -- 'gurobi': uses Gurobi as the solver, for more information
%                          see "http://www.gurobi.com"
% Output:
%   Decision -- the output is one of the following messages:
%               - 'Input bounds are not satisfied'
%               - 'The model is invalidated'
%               - 'The model is not invalidated'
%   sol -- A solution structure that is the output of YALMIP
%
% Syntax:
%   Decision = SWA_MILP(SYS,y,u,p,u_bound,x_bound,mn_bound, solver);
%
% Author: MI4Hybrid
% Date: July 13th, 2015


%% Input, States, Output dimensions and time horizon 
[n_y N] = size(y);
n_u = size(u);
num_m = size(SYS.mode,2);   % number of modes
n = size(SYS.mode(1).A,1);

%% Initiate system modes
for i = 1: num_m
    Mode(i).A = SYS.mode(i).A;
    Mode(i).B = SYS.mode(i).B;
    Mode(i).C = SYS.mode(i).C;
    Mode(i).D = SYS.mode(i).D;
    Mode(i).f = SYS.mode(i).f;
    Mode(i).g = SYS.mode(i).g;
end

M = x_bound;
eps = mn_bound;
U = u_bound;


%% Checking input bounds
for i = 1:n_u
    u_norm(1,i) = norm(u(i,:),p);
end
if sum(u_norm > u_bound) >0
    Decision = 'Input bounds are not satisfied';
    return
end


%% Building the constraints
for i = 1:num_m
    sys(i).g = [-Mode(i).C -ones(n_y,1)];
    sys(i).q = [-Mode(i).D ones(n_y,1)];
    sys(i).h = [-Mode(i).A eye(n)];
    sys(i).l = [-Mode(i).B];
end

s = binvar(N,num_m);
e = sdpvar(num_m*N,n_y);
X = sdpvar(num_m*N,2*n);

Constraint = [];
for t = 1:N-1    % time index
    
    for sys_ind = 1: num_m  % mode index
        
        for j = 1:n_y   % constraints for output equation
            Constraint = [Constraint sys(sys_ind).g(j,:)*[X((t-1)* ...
                num_m+sys_ind,1:n) e((t-1)*num_m+sys_ind,j)]'+ ...
                s(t,sys_ind)*sys(sys_ind).q(j,:)*[u(:,t) y(j,t)]' ==0];
        end
        
        for state_ind = 1:n  % constraints for state equation
            Constraint = [Constraint sys(sys_ind).h(state_ind,:)* ...
                X((t-1)*num_m+sys_ind,:)'+ ...
                s(t,sys_ind)*sys(sys_ind).l(state_ind,:)*u(:,t)'==0];
            Constraint = [Constraint norm(X((t-1)*num_m+sys_ind, ...
                state_ind),inf)<=M(state_ind)*s(t,sys_ind)];
        end
        
        
        Constraint = [Constraint norm(e((t-1)*num_m+sys_ind,:),inf)<=...
            eps*s(t,sys_ind)];
    end
    Constraint = [Constraint sum(s(t,:))==1];
    Constraint = [Constraint sum(X((t-1)*num_m+1:t*num_m,n+1:2*n),1)...
        ==sum(X(t*num_m+1:(t+1)*num_m,1:n),1)];
end

% Constraint for the last output equation
for sys_ind = 1: num_m
    for j = 1:n_y   % constraints for output equation
        Constraint = [Constraint sys(sys_ind).g(j,:)*[X((N-1)*num_m+ ...
            sys_ind,n+1:2*n) e((N-1)*num_m+sys_ind,j)]'+ ...
            s(N,sys_ind)*sys(sys_ind).q(j,:)*[u(:,N) y(j,N)]' ==0];
    end
    Constraint = [Constraint norm(e((N-1)*num_m+sys_ind,:),inf)<= ...
        eps*s(N,sys_ind)];
    Constraint = [Constraint norm(X((N-1)*num_m+sys_ind,n+1:2*n),inf)<=...
        M(sys_ind)*s(N,sys_ind)];
end
Constraint = [Constraint sum(s(N,:))==1];


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



