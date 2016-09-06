function [Decision,sol,z_c,z_f] = Tdetect_swa_milp(SYS,SYS_f, ...
    T,mn_bnd,pn_bnd,input_low, input_high, state_low, state_high,solver)

% This function implements MILP-based T-detectability approach for the
% switched state-space models subject to measurement and process noise, but
% without uncertainty in parameters discussed in "Guaranteed Model-Based Fault 
% Detection in Cyber-Physical Systems: A Model Invalidation Approach" by
% F. Harirchi and N. Ozay.
%
% Preleminaries for this function are: YALMIP and (cplex or gurobi).
%
% Inputs/Parameters:
%   SYS -- a user-defined system class (see StateSpace.m)
%   SYS_f -- a user-defined fault class (see StateSpace.m)
%   mn_bound -- an n_y-by-2 matrix specifying the bound for measurement
%               noise for system (1st column) and fault (2nd column) where
%               n_y is the dimension of the output.
%   pn_bound -- an n-by-2 matrix specifying the bound for process
%               noise for the system (1st column) and fault (2nd column) where
%               n is the dimension of the states.
%   input_high -- an n_u-by-1 vector specifying upperbound on each element
%                  of input, where n_u is the number of inputs.
%   input_low -- an n_u-by-1 vector specifying lowerbound on each element
%                  of input, where n_u is the number of inputs.
%   state_high -- an n-by-2 matrix specifying upperbound on each element of 
%                  state variables for the system (1st column) and fault 
%                  (2nd column) where n is the number of states.
%   state_low -- an n-by-2 matrix specifying lowerbound on each element of 
%                  state variables for the system (1st column) and fault 
%                  (2nd column) where n is the number of states.
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
% Date: January 13th, 2016

%% Check if the system model is valid for this function
if(strcmp(SYS.mark,'ss')~=1&&strcmp(SYS.mark,'swss')~=1)
    error('The system model must be a state-space model.');
end

%% Check if the fault model is valid for this function
if(strcmp(SYS_f.mark,'ss')~=1&&strcmp(SYS_f.mark,'swss')~=1)
    error('The fault model must be a state-space model.');
end


%% Input, States, Output dimensions and time horizon
n_y = size(SYS.mode(1).C,1); % output dimension and time horizon
n_u = size(SYS.mode(1).B,2); % input dimension
n_mode = size(SYS.mode,2); % number of system modes
nf_mode = size(SYS_f.mode,2); % number of fault modes
n = size(SYS.mode(1).A,1); % state dimension


%% Set up default values for empty paramters
if(isempty(mn_bnd))
    mn_bnd=zeros(n_y,2);
end
if(isempty(pn_bnd))
    pn_bnd=zeros(n,2);
end
if(isempty(input_low))
    input_low=zeros(n_u,1)-inf;
end
if(isempty(input_high))
    input_high=zeros(n_u,1)+inf;
end
if(isempty(state_low))
    state_low=zeros(n,2)-inf*ones(n,2);
end
if(isempty(state_high))
    state_high=zeros(n,2)+inf*ones(n,2);
end
if(isempty(solver))
    solver='cplex';
end

%% Convert scalars to vectors
if(size(mn_bnd,1)==1&&n_y>1)
    if(size(mn_bnd,2)==1)
        mn_bnd=ones(n_y,2)*mn_bnd;
            warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
    elseif (size(mn_bnd,2)==2)
        mn_bnd=[ones(n_y,1)*mn_bnd(1,1) ones(n_y,1)*mn_bnd(1,2)];
            warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
    else
    error(['mn_bnd dimensions are not consistent.']);
    end
end

if(size(pn_bnd,1)==1&&n>1)
    if(size(pn_bnd,2)==1)
        pn_bnd=ones(n,2)*pn_bnd;
            warning(['Bound for process noise is a scalar, converted to'...
        ' a vector with identical entries.']);
    elseif (size(pn_bnd,2)==2)
        pn_bnd=[ones(n,1)*pn_bnd(1,1) ones(n,1)*pn_bnd(1,2)];
            warning(['Bound for process noise is a scalar, converted to'...
        ' a vector with identical entries.']);
    else
    error(['pn_bnd dimensions are not consistent.']);
    end
end

if(size(input_low,1)==1&&n_u>1)
    if(size(input_low,2)==1)
        input_low=ones(n_u,2)*input_low;
            warning(['Lowerbound for input is a scalar, converted to'...
        ' a vector with identical entries.']);
    elseif (size(input_low,2)==2)
        input_low=[ones(n_u,1)*input_low(1,1) ones(n_u,1)*input_low(1,2)];
            warning(['Lowerbound for input is a scalar, converted to'...
        ' a vector with identical entries.']);
    else
    error(['input_low dimensions are not consistent.']);
    end
end

if(size(input_high,1)==1&&n_u>1)
    if(size(input_high,2)==1)
        input_high=ones(n_u,2)*input_high;
            warning(['Upperbound for input is a scalar, converted to'...
        ' a vector with identical entries.']);
    elseif (size(input_high,2)==2)
        input_high=[ones(n_u,1)*input_high(1,1) ones(n_u,1)*input_high(1,2)];
            warning(['Upperbound for input is a scalar, converted to'...
        ' a vector with identical entries.']);
    else
    error(['input_high dimensions are not consistent.']);
    end
end


if(size(state_low,1)==1&&n>1)
    if(size(state_low,2)==1)
        state_low=ones(n,2)*state_low;
        warning(['State lower bound is a scalar, converted to a vector with '...
        'identical entries.']);
    elseif (size(state_low,2)==2)
        state_low=[ones(n,1)*state_low(1,1) ones(n,1)*state_low(1,2)];
        warning(['State lower bound is a scalar, converted to a vector with '...
        'identical entries.']);
    else
    error(['state_low dimensions are not consistent.']);
    end
end

if(size(state_high,1)==1&&n>1)
    if(size(state_high,2)==1)
        state_high=ones(n,2)*state_high;
        warning(['State upperbound is a scalar, converted to a vector with '...
        'identical entries.']);
    elseif (size(state_high,2)==2)
        state_high=[ones(n,1)*state_high(1,1) ones(n,1)*state_high(1,2)];
        warning(['State upperbound is a scalar, converted to a vector with '...
        'identical entries.']);
    else
    error(['state_up dimensions are not consistent.']);
    end
end

%% Check the bounds
if(size(mn_bnd,1)~=n_y||size(mn_bnd,2)~=2)
    error('The number of bounds for measurement noise is not correct.');
end
if(size(pn_bnd,1)~=n||size(pn_bnd,2)~=2)
    error('The number of bounds for process noise is not correct.');
end
if(size(input_low,1)~=n_u||~isvector(input_low))
    error('The number of lowerbounds for inputs is not correct.');
end
if(size(input_high,1)~=n_u||~isvector(input_high))
    error('The number of upperbounds for inputs is not correct.');
end
if(size(state_low,1)~=n||size(state_low,2)~=2)
    error('The number of lowerbounds for states is not correct.');
end
if(size(state_high,1)~=n||size(state_high,2)~=2)
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

%% Initiate fault modes
for i = 1: nf_mode
    fMode(i).A = SYS_f.mode(i).A;
    fMode(i).B = SYS_f.mode(i).B;
    fMode(i).C = SYS_f.mode(i).C;
    fMode(i).D = SYS_f.mode(i).D;
    fMode(i).f = SYS_f.mode(i).f;
    fMode(i).g = SYS_f.mode(i).g;
end

M_u = state_high(:,1);
M_l = state_low(:,1);
Mf_u = state_high(:,2);
Mf_l = state_low(:,2);
eps = mn_bnd(:,1);
epsf = mn_bnd(:,2);
PNB = pn_bnd(:,1);
PNBf = pn_bnd(:,2);
U_u = input_high(:,1);
U_l = input_low(:,1);

%% Optimization variables
d = binvar(n_mode,nf_mode,T);      % binary variables
z_c = sdpvar(n,n_mode,nf_mode,T);  % current state
z_f = sdpvar(n,n_mode,nf_mode,T);  % future state
zf_c = sdpvar(n,n_mode,nf_mode,T);  % fault current state
zf_f = sdpvar(n,n_mode,nf_mode,T);  % fault future state
theta = sdpvar(n_y,n_mode,nf_mode,T); % measurement noise (system)
thetaf = sdpvar(n_y,n_mode,nf_mode,T); % measurement noise (fault)
eta = sdpvar(n,n_mode,nf_mode,T); % process noise (system)
etaf = sdpvar(n,n_mode,nf_mode,T); % process noise (fault)
tt = sdpvar(n_u,n_mode,nf_mode,T); % input 

Constraint = [];

% state and output equation constraints:
for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for flt_ind = 1: nf_mode  % fault mode index
            
         Constraint = [Constraint z_f(:,sys_ind,flt_ind,t) - ...
             Mode(sys_ind).A * z_c(:,sys_ind,flt_ind,t) - Mode(sys_ind).B...
             * tt(:,sys_ind,flt_ind,t) - eta(:,sys_ind,flt_ind,t)   == ...
             d(sys_ind,flt_ind,t)* Mode(sys_ind).f];
             
         Constraint = [Constraint zf_f(:,sys_ind,flt_ind,t) - ...
             fMode(flt_ind).A * zf_c(:,sys_ind,flt_ind,t) - fMode(flt_ind).B...
             * tt(:,sys_ind,flt_ind,t)- etaf(:,sys_ind,flt_ind,t) == ...
             d(sys_ind,flt_ind,t)* (fMode(flt_ind).f)];
            
         Constraint = [Constraint Mode(sys_ind).C * z_f(:,sys_ind,flt_ind,t)...
             + theta(:,sys_ind,flt_ind,t)+ d(sys_ind,flt_ind,t)*Mode(sys_ind).g ...
             ==  fMode(flt_ind).C * zf_f(:,sys_ind,flt_ind,t) + ...
             thetaf(:,sys_ind,flt_ind,t)+ d(sys_ind,flt_ind,t)*fMode(flt_ind).g];
        end
    end
end

% integer constraints
for t = 1:T
    Constraint = [Constraint sum(sum(d(:,:,t)))==1];
end

% same initial condition constraint
Constraint = [Constraint sum(sum(z_c(:,:,:,1),3),2)==...
    sum(sum(zf_c(:,:,:,1),3),2)];

% current and future state constraints
for t = 1:T-1
    for flt_ind = 1:nf_mode
        Constraint = [Constraint sum(z_c(:,:,flt_ind,t+1),2)==...
            sum(z_f(:,:,flt_ind,t),2)];
    end

    for sys_ind = 1:n_mode
        Constraint = [Constraint sum(zf_c(:,sys_ind,:,t+1),3)==...
            sum(zf_f(:,sys_ind,:,t),3)];
    end
end
% admissible set constraints
for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for flt_ind = 1: nf_mode  % fault mode index
            for i = 1:n_y
            Constraint = [Constraint  norm(theta(i,sys_ind,flt_ind,t),inf)<=...
                d(sys_ind,flt_ind,t)* eps(i)];
            Constraint = [Constraint  norm(thetaf(i,sys_ind,flt_ind,t),inf)<=...
                d(sys_ind,flt_ind,t)* epsf(i)];
            end
            for i = 1:n
            Constraint = [Constraint  norm(eta(i,sys_ind,flt_ind,t),inf)<=...
                d(sys_ind,flt_ind,t)* PNB(i)];
            Constraint = [Constraint  norm(etaf(i,sys_ind,flt_ind,t),inf)<=...
                d(sys_ind,flt_ind,t)* PNBf(i)];
            end
        end
    end
end

for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for flt_ind = 1: nf_mode  % fault mode index
        for i = 1:n
            Constraint = [Constraint d(sys_ind,flt_ind,t)* (M_l(i))<=...
                z_c(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint z_c(i,sys_ind,flt_ind,t)<= ...
                d(sys_ind,flt_ind,t)*M_u(i)];
            Constraint = [Constraint d(sys_ind,flt_ind,t)* (Mf_l(i))<=...
                zf_c(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint zf_c(i,sys_ind,flt_ind,t)<=...
                d(sys_ind,flt_ind,t)*Mf_u(i)];
        end
        
        for i = 1:n_u
            Constraint = [Constraint  d(sys_ind,flt_ind,t)* (U_l(i))<=...
                tt(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint tt(i,sys_ind,flt_ind,t)<= ...
                d(sys_ind,flt_ind,t)*U_u(i)];
        end
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
        Decision = ['The fault is ', num2str(T), '-detectable for the system'];
    else
        Decision = ['The fault is not ', num2str(T), '-detectable for the system'];
    end
end
if strcmp(solver, 'gurobi')
    if strcmp(sol.info, 'Infeasible problem (GUROBI-GUROBI)')
        Decision = ['The fault is ', num2str(T), '-detectable for the system'];
    else
        Decision = ['The fault is not ', num2str(T), '-detectable for the system'];
    end
end