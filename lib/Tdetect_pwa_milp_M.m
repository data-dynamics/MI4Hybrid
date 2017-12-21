function [Decision,sol] = Tdetect_pwa_milp_M(SYS,SYS_f, ...
    T,mn_bound,pn_bound,input_bound, state_bound,solver)

% This function implements MILP-based model invalidation approach for the
% piecewise switch affine(PWA) system.
%
% Preleminaries for this function are: YALMIP and (cplex or gurobi).
%
% Inputs/Parameters:
%   SYS -- a user-defined system class (see StateSpace.m)
%   SYS_f -- a user-defined system class (see StateSpace.m)
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
%               - 'The fault is T-detectable'
%               - 'The fault is not T-detectable'
%   sol -- A solution structure that is the output of YALMIP
%
% Syntax:
%   [Decision,sol] = Tdetect_swa_milp(SYS,SYS_f, ...
%    T,mn_bnd,pn_bnd,input_low, input_high, state_low, state_high,solver)
%
% Author: MI4Hybrid
% Modified: Supratim Ghosh
% Date: December 21st, 2017

%% Check if the system model is valid for this function
if(strcmp(SYS.mark,'pwa')~=1)
    error('he system model must be a piecewise switch affine(PWA) model..');
end

%% Check if the fault model is valid for this function
if(strcmp(SYS_f.mark,'pwa')~=1)
    error('he system model must be a piecewise switch affine(PWA) model.');
end

%% Input, States, Output dimensions and time horizon
n_y = size(SYS.mode(1).C,1); % output dimension and time horizon
n_u = size(SYS.mode(1).B,2); % input dimension
n_mode = size(SYS.mode,2); % number of system modes
nf_mode = size(SYS_f.mode,2); % number of fault modes
n = size(SYS.mode(1).A,1); % state dimension
d = size(SYS.mode(1).P,1); % dimension of hyper-plane constraints
%% Set up default values for empty paramters
if(isempty(pn_bound))
    pn_bound=zeros(n,2);
end
if(isempty(mn_bound))
    mn_bound=zeros(n_y,2);
end
if(isempty(input_bound))
    input_bound=zeros(n_i,1)+inf;
end
if(isempty(state_bound))
    state_bound=zeros(n,2)+inf;
end
if(isempty(solver))
    solver='gurobi';
end


%% Convert scalars to vectors

%% Convert measurement noise bound to a vector
if(size(mn_bound,1)==1&&n_y>1)
    if(size(mn_bound,2)==1)
        mn_bound=ones(n_y,2)*mn_bound;
            warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
    elseif (size(mn_bound,2)==2)
        mn_bound=[ones(n_y,1)*mn_bound(1,1) ones(n_y,1)*mn_bound(1,2)];
            warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
    else
    error(['mn_bound dimensions are not consistent.']);
    end
elseif (size(mn_bound,2)==1)
    mn_bound=ones(n_y,2)*mn_bound;
end

%% Convert process noise bound to a vector

if(size(pn_bound,1)==1&&n>1)
    if(size(pn_bound,2)==1)
        pn_bound=ones(n,2)*pn_bound;
            warning(['Bound for process noise is a scalar, converted to'...
        ' a vector with identical entries.']);
    elseif (size(pn_bound,2)==2)
        pn_bound=[ones(n,1)*pn_bound(1,1) ones(n,1)*pn_bound(1,2)];
            warning(['Bound for process noise is a scalar, converted to'...
        ' a vector with identical entries.']);
    else
    error(['pn_bound dimensions are not consistent.']);
    end
end

%% Convert input bound to a vector

if(size(input_bound,1)==1&&n_u>1)
    if(size(input_bound,2)==1)
        input_bound=ones(n_u,2)*input_bound;
            warning(['Lowerbound for input is a scalar, converted to'...
        ' a vector with identical entries.']);
    elseif (size(input_bound,2)==2)
        input_bound=[ones(n_u,1)*input_bound(1,1) ones(n_u,1)*input_bound(1,2)];
            warning(['Lowerbound for input is a scalar, converted to'...
        ' a vector with identical entries.']);
    else
    error(['input_low dimensions are not consistent.']);
    end
end

%% Convert state bound to a vector

if(size(state_bound,1)==1&&n>1)
    if(size(state_bound,2)==1)
        state_bound=ones(n,2)*state_bound;
        warning(['State lower bound is a scalar, converted to a vector with '...
        'identical entries.']);
    elseif (size(state_bound,2)==2)
        state_bound=[ones(n,1)*state_bound(1,1) ones(n,1)*state_bound(1,2)];
        warning(['State lower bound is a scalar, converted to a vector with '...
        'identical entries.']);
    else
    error(['state_low dimensions are not consistent.']);
    end
end


%% Check the bounds
if(size(mn_bound,1)~=n_y||size(mn_bound,2)~=2)
    error('The number of bounds for measurement noise is not correct.');
end
if(size(pn_bound,1)~=n||size(pn_bound,2)~=2)
    error('The number of bounds for process noise is not correct.');
end
% if(size(input_bound,1)~=n_u||~isvector(input_bound))
%     error('The number of bounds for inputs is not correct.');
% end
% if(size(input_high,1)~=n_u||~isvector(input_high))
%     error('The number of upperbounds for inputs is not correct.');
% end
if(size(state_bound,1)~=n||size(state_bound,2)~=2)
    error('The number of lowerbounds for states is not correct.');
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

%% Initiate fault modes
for i = 1: nf_mode
    fMode(i).A = SYS_f.mode(i).A;
    fMode(i).B = SYS_f.mode(i).B;
    fMode(i).C = SYS_f.mode(i).C;
    fMode(i).D = SYS_f.mode(i).D;
    fMode(i).f = SYS_f.mode(i).f;
    fMode(i).g = SYS_f.mode(i).g;
    fMode(i).P = SYS_f.mode(i).P;
    fMode(i).M = SYS_f.mode(i).M;
end

% Specify the bounds of various quantities
S = state_bound(:,1);
Sf = state_bound(:,2);
eps = mn_bound(:,1);
epsf = mn_bound(:,2);
gamma = pn_bound(:,1);
gammaf = pn_bound(:,2);
U = input_bound;

%% Building the constraints


a = binvar(n_mode,nf_mode,T,'full'); %sequence matrix

% real variables
x = sdpvar(n,T,'full');  % current state
xf = sdpvar(n,T,'full');  % fault state
eta = sdpvar(n_y,T,'full'); % measurement noise (system)
etaf = sdpvar(n_y,T,'full'); % measurement noise (fault)
nu = sdpvar(n,T,'full'); % process noise (system)
nuf = sdpvar(n,T,'full'); % process noise (fault)
u = sdpvar(n_u,T,'full'); % input 
sdpvar delta

%% Creating Constraints
Constraint = [];
M1 = 1000*ones(d,1); % big M bound for hyperplane constraints
M2 = 10000*ones(n,1); % big M bound for state constraints
M3 = 10000*ones(n_y,1); % big M bound for output constraints

for t = 1:T-1    % time index
    for sys_ind = 1: n_mode  % system mode index 
        for flt_ind = 1: nf_mode  % system mode index
            for i = 1:n
                % State Constraints (system)
                Constraint = [Constraint, -M2(i)*(1-a(sys_ind,flt_ind,t)) <= x(i,t+1) - Mode(sys_ind).A(i,:) * x(:,t)-Mode(sys_ind).B(i,:)*u(:,t)-nu(i,t)-Mode(sys_ind).f(i)];
                Constraint = [Constraint, M2(i)*(1-a(sys_ind,flt_ind,t)) >= x(i,t+1) - Mode(sys_ind).A(i,:) * x(:,t)-Mode(sys_ind).B(i,:)*u(:,t)-nu(i,t)-Mode(sys_ind).f(i)];
            end
            for i = 1:n
                % State Constraints (fault)
                Constraint = [Constraint, -M2(i)*(1-a(sys_ind,flt_ind,t)) <= xf(i,t+1)-fMode(flt_ind).A(i,:)*xf(:,t)-fMode(flt_ind).B(i,:)*u(:,t)-nuf(i,t)-fMode(flt_ind).f(i)];
                Constraint = [Constraint, M2(i)*(1-a(sys_ind,flt_ind,t)) >= xf(i,t+1)-fMode(flt_ind).A(i,:)*xf(:,t)-fMode(flt_ind).B(i,:)*u(:,t)-nuf(i,t)-fMode(flt_ind).f(i)];
            end
            for i = 1:n_y
                % Output Constraints
                Constraint = [Constraint, -M3(i)*(1-a(sys_ind,flt_ind,t)) <= Mode(sys_ind).C(i,:)*x(:,t)+Mode(sys_ind).D(i,:)*u(:,t)+eta(i,t)+Mode(sys_ind).g(i)-fMode(flt_ind).C(i,:)*xf(:,t)-fMode(flt_ind).D(i,:)*u(:,t)-etaf(i,t)-fMode(flt_ind).g(i)];
                Constraint = [Constraint, M3(i)*(1-a(sys_ind,flt_ind,t)) >= Mode(sys_ind).C(i,:)*x(:,t)+Mode(sys_ind).D(i,:)*u(:,t)+eta(i,t)+Mode(sys_ind).g(i)-fMode(flt_ind).C(i,:)*xf(:,t)-fMode(flt_ind).D(i,:)*u(:,t)-etaf(i,t)-fMode(flt_ind).g(i)];
            end
        end
    end
end

%% integer constraints
for t = 1:T
    Constraint = [Constraint sum(sum(a(:,:,t)))==1];
end
             
for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for flt_ind = 1: nf_mode  % fault mode index
            Constraint = [Constraint, Mode(sys_ind).P*x(:,t)+Mode(sys_ind).M<= M1*(1-a(sys_ind,flt_ind,t))];
            Constraint = [Constraint, fMode(flt_ind).P*xf(:,t)+fMode(flt_ind).M<= M1*(1-a(sys_ind,flt_ind,t))];
        end
    end
end
                   
%% admissible set constraints
% Satisfaction of noise bounds
for t = 1:T    % time index
        for i = 1:n_y
            Constraint = [Constraint  norm(eta(i,t),inf)<=eps(i)];
            Constraint = [Constraint,  norm(etaf(i,t),inf)<=epsf(i)];
        end
        for i = 1:n
            Constraint = [Constraint  norm(nu(i,t),inf)<=gamma(i)];
            Constraint = [Constraint,  norm(nuf(i,t),inf) <= gammaf(i)];
        end
end

% Satisfaction of input and state bounds
for t = 1:T    % time index
    for i = 1:n
        Constraint = [Constraint, (-S(i))<= x(i,t)];
        Constraint = [Constraint, x(i,t)<= S(i)];
        Constraint = [Constraint, (-Sf(i)) <= xf(i,t)];
        Constraint = [Constraint, xf(i,t) <= Sf(i)];
    end
    
    for i = 1:n_u
        Constraint = [Constraint, (-U(i))<= u(i,t)];
        Constraint = [Constraint, u(i,t) <= U(i)];
    end
end

% obj = delta;
obj = [];

%% Setting up Options
options = sdpsettings('verbose',1,'solver',solver);
options.gurobi.Heuristics=0.7;
%% Solve the problem
sol = optimize(Constraint,obj,options);
sol.info
value(sol)
% AD = value(a)
% UD = value(u)
% XD = value(x)
% XDf = value(xf)
%% Print Results
if strcmp(solver,'cplex')
    if strcmp (sol.info , 'Infeasible problem (CPLEX-IBM)')
        Decision = ['The fault is ', num2str(T), '-detectable for the system'];
    else
        Decision = ['The fault is not ', num2str(T), '-detectable for the system'];
    end
end
if strcmp(solver, 'gurobi')
    if (strcmp(sol.info, 'Infeasible problem (GUROBI-GUROBI)')|| strcmp(sol.info,'Either infeasible or unbounded (GUROBI-GUROBI)'))
        Decision = ['The fault is ', num2str(T), '-detectable for the system'];
    else
        Decision = ['The fault is not ', num2str(T), '-detectable for the system'];
    end
end

