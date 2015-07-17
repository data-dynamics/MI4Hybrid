function [Decision, sol,Constraint,e,s] = SARX_MILP(SYS,y,u,p,u_bnd,mn_bnd,solver)

% This function implements MILP-based model invalidation approach for 
% switched state-space models discussed in "Model Invalidation for Switched 
% Affine Systems with Applications to Fault and Anomaly Detection" in 
% ADHS15 by F. Harirchi and N. Ozay.
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
%   u_bnd -- an n_u-by-1 vector specifying p-norm bound on input values, 
%               where n_u is number of inputs.
%   mn_bnd -- an n_y-by-1 vector specifying the bound for process noise
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
%   Decision = SARX_MILP(SYS,y,u,p,u_bnd,mn_bnd, solver);
%
% Author: MI4Hybrid
% Date: July 15th, 2015


%% Input, States, Output dimensions and time horizon 
[n_y N] = size(y);          % number of outputs and time horizon
n_u = size(u);    % number of inputs
% num_m = 1;
num_m = size(SYS.mode,2);   % number of modes

%% Initiate system modes
for i = 1: num_m
    Mode(i).A = SYS.mode(i).A;
    Mode(i).C = SYS.mode(i).C;
    Mode(i).f = SYS.mode(i).f;
end

eps = mn_bnd;
U = u_bnd;

%% Checking input bounds
for i = 1:n_u
    u_norm(1,i) = norm(u(i,:),p);
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
for i = 1:num_m
    MODE(i).A(:,:,1) = -eye(n_y);
    for k = 1:n_a
        MODE(i).A(:,:,k+1) = Mode(i).A(:,:,k);
    end
end

for i = 1:num_m
    sys(i).h = [];
    for j = 1: n_y
        for k = 1:n_a+1
            sys(i).h = [sys(i).h MODE(i).A(:,j,k)];
        end
    end
    sys(i).g = [];
    for j = 1:n_u
        for k = 1:n_c
            sys(i).g = [sys(i).g Mode(i).C(:,j,k)];
        end
    end
    sys(i).f = Mode(i).f;
end


% Define variables
s = binvar(N,num_m);
e = sdpvar(n_y*(n_a+1),num_m*N);


Constraint = [];
for t = n+1:N-1   % time index
    
    for sys_ind = 1: num_m  % mode index
        
        for j = 1:n_y   % constraints for output equation
           Constraint = [Constraint -sys(sys_ind).h(j,:)*e(:,(t-1)*num_m...
           +sys_ind)+s(t,sys_ind)*(sys(sys_ind).g(j,:)*reshape(u(:,...
           t-1:-1:t-n_c)',1,[])'+sys(sys_ind).h(j,:)*reshape(y(:,...
           t:-1:t-n_a)',1,[])')+ sys(sys_ind).f(j) ==0];
        end
        
        for i=1:n_y
            %constraints for noise variables
            Constraint = [Constraint norm(e((i-1)*(n_a+1)+1:i*(n_a+1),(t-1)*num_m+sys_ind),inf)<=...
             eps(n_y)*s(t,sys_ind)];
%             Constraint = [Constraint norm(e(:,(t-1)*num_m+sys_ind),inf)<=...
%              eps(1)*s(t,sys_ind)];
        end
                
    end
    Constraint = [Constraint sum(s(t,:))==1];
                
    % constraints for noise variables
    for i=1:n_y
            Constraint = [Constraint sum(e((i-1)*(n_a+1)+2:(i)*(n_a+1),(t-1)*num_m+1:...
            t*num_m),2)==sum(e((i-1)*(n_a+1)+1:i*(n_a+1)-1,(t-2)*num_m+1:(t-1)*num_m),2)];   
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



