function [Decision,sol,z_c,z_f] = wTdetect_uswa_milp(SYS,SYS_f,Delta,Delta_f, ...
    T,W,mn_bound,pn_bound,input_bound,state_bound, un_bound, solver)

% This function implements MILP-based model invalidation approach for the
% switched state-space models discussed in "Model Invalidation for Switched
% Affine Systems with Applications to Fault and Anomaly Detection" by
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
%   input_bound -- an n_u-by-2 matrix specifying upperbound on absolute value
%                  of input for the system (1st column) and fault (2nd column)
%                  where n_u is number of inputs
%   state_bound -- an n-by-2 vector specifying infinity norm bound on state
%                  variables for the system (1st column) and fault (2nd column)
%                  where n is number of states
%   uncertainty_bound -- 1-by-2 vector, each value corresponds to the bound
%                        on uncertainty
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
%   Decision = Tdetect_swa_milp(SYS,SYS_f,mn_bound, pn_bound,input_bound,state_bound,solver);
%
% Author: MI4Hybrid
% Date: January 13th, 2016

%% Check if the system model is valid for this function
if(strcmp(SYS.mark,'ss')~=1&&strcmp(SYS.mark,'swss')~=1)
    error('The system model must be a state-space model.');
end

%% Check if the fault model is valid for this function
if(strcmp(SYS_f.mark,'ss')~=1&&strcmp(SYS_f.mark,'swss')~=1)
    error('The system model must be a state-space model.');
end

%% Check if the input, output, and the model are consistent
% if(length(input)~=length(output))
%     error('The input length and output length are not the same.');
% end
% if(size(SYS_f.mode(1).B,2)~=size(SYS.mode(1).B,2))
%     error('The input of the two models is not consistent.');
% end
% if(size(SYS.mode(1).C,1)~= size(SYS_f.mode(1).C,1))
%     error('The output of the two models is not consistent.');
% end
% if(size(SYS.mode(1).A,1)~= size(SYS_f.mode(1).A,1))
%     error('The states of the two models is not consistent.');
% end

%% Input, States, Output dimensions and time horizon
n_y = size(SYS.mode(1).C,1); % output dimension and time horizon
n_u = size(SYS.mode(1).B,2); % input dimension
n_mode = size(SYS.mode,2); % number of modes
nf_mode = size(SYS_f.mode,2); % number of modes
n = size(SYS.mode(1).A,1); % state dimension

%% Use the default solver if it is not specified
if(nargin==6)
    solver='cplex';
end

% %% Set up default values for empty paramters
% if(isempty(mn_bound))
%     mn_bound=zeros(n_y,1);
% end
% if(isempty(input_bound))
%     input_bound=zeros(n_i,1)+inf;
% end
% if(isempty(state_bound))
%     state_bound=zeros(n,1)+inf;
% end
if(isempty(solver))
    solver='cplex';
end

%% Convert scalars to vectors
% if(length(mn_bound)==1&&n_y>1)
%     mn_bound=ones(n_y,1)*mn_bound;
%     warning(['Bound for measurement noise is a scalar, converted to'...
%         ' a vector with identical entries.']);
% end
% if(length(input_bound)==1&&n_i>1)
%     input_bound=ones(n_i,1)*input_bound;
%     warning(['Input bound is a scalar, converted to a vector with '...
%         'identical entries.']);
% end
% if(length(state_bound)==1&&n>1)
%     state_bound=ones(n,1)*state_bound;
%     warning(['State bound is a scalar, converted to a vector with '...
%         'identical entries.']);
% end

% %% Check the bounds
% if(length(mn_bound)~=n_y||~isvector(mn_bound))
%     error('The number of bounds for measurement noise is not correct.');
% end
% if(length(input_bound)~=n_i||~isvector(input_bound))
%     error('The number of bounds for inputs is not correct.');
% end
% if(length(state_bound)~=n||~isvector(state_bound))
%     error('The number of bounds for states is not correct.');
% end

%% Initiate system modes
for i = 1: n_mode
    Mode(i).A = SYS.mode(i).A;
    Mode(i).B = SYS.mode(i).B;
    Mode(i).C = SYS.mode(i).C;
    Mode(i).D = SYS.mode(i).D;
    Mode(i).f = SYS.mode(i).f;
    Mode(i).g = SYS.mode(i).g;
end

%% Initiate system modes
for i = 1: nf_mode
    fMode(i).A = SYS_f.mode(i).A;
    fMode(i).B = SYS_f.mode(i).B;
    fMode(i).C = SYS_f.mode(i).C;
    fMode(i).D = SYS_f.mode(i).D;
    fMode(i).f = SYS_f.mode(i).f;
    fMode(i).g = SYS_f.mode(i).g;
end

M_u = state_bound(1,1);
M_l = state_bound(2,1);
Mf_u = state_bound(1,2);
Mf_l = state_bound(2,2);
eps = mn_bound(:,1);
epsf = mn_bound(:,2);
U_u = input_bound(1,1);
U_l = input_bound(2,1);
Uf_u = input_bound(1,2);
Uf_l = input_bound(2,2);
P = pn_bound(:,1);
Pf = pn_bound(:,2);
delta = un_bound(:,1);
deltaf = un_bound(:,2);

% %% Checking input bounds
% for i = 1:n_i
%     u_norm(1,i) = norm(input(i,:),SYS.input_norm(i));
% end
% if sum(u_norm > input_bound) >0
%     Decision = 'Input bounds are not satisfied';
%     return
% end

%% Building the constraints
% for i = 1:n_mode
%     sys(i).g = [-Mode(i).C -ones(n_y,1)];
%     sys(i).q = [-Mode(i).D ones(n_y,1)];
%     sys(i).h = [-Mode(i).A eye(n)];
%     sys(i).l = [-Mode(i).B];
% end

d = binvar(n_mode,nf_mode,T);
z_c = sdpvar(n,n_mode,nf_mode,T);  % current
z_f = sdpvar(n,n_mode,nf_mode,T);  % future
zf_c = sdpvar(n,n_mode,nf_mode,T);  % current
zf_f = sdpvar(n,n_mode,nf_mode,T);  % future
theta = sdpvar(n_y,n_mode,nf_mode,T);
thetaf = sdpvar(n_y,n_mode,nf_mode,T);
tt = sdpvar(n_u,n_mode,nf_mode,T);
Deltaf = sdpvar(n_mode,nf_mode,T);
Deltaf_f = sdpvar(n_mode,nf_mode,T);
% Delta_Ax = sdpvar(n,n,n_mode,nf_mode,T);
% Deltaf_Ax = sdpvar(n,n,n_mode,nf_mode,T);
% Delta_Cx = sdpvar(n_y,n,n_mode,nf_mode,T);
% Deltaf_Cx = sdpvar(n_y,n,n_mode,nf_mode,T);
% Delta_Bu = sdpvar(n,n_u,n_mode,nf_mode,T);
% Deltaf_Bu = sdpvar(n,n_u,n_mode,nf_mode,T);



Constraint = [];


% for t = 1:T-1    % time index
%     for sys_ind = 1: n_mode  % system mode index
%         for flt_ind = 1: nf_mode  % fault mode index
%              Constraint = [Constraint, z_f(:,sys_ind,flt_ind,t) - Mode(sys_ind).A * z_c(:,sys_ind,flt_ind,t)  ...
%                 - sum(Delta_Ax(:,:,sys_ind,flt_ind,t),2)- Mode(sys_ind).B * tt(:,sys_ind,flt_ind,t) - ...
%                  sum(Delta_Bu(:,:,sys_ind,flt_ind,t),2)- d(sys_ind,flt_ind,t)* Mode(sys_ind).g == zeros(n,1)];
%              Constraint = [Constraint, zf_f(:,sys_ind,flt_ind,t) - fMode(flt_ind).A * zf_c(:,sys_ind,flt_ind,t)  ...
%                 - sum(Deltaf_Ax(:,:,sys_ind,flt_ind,t),2)- fMode(flt_ind).B * tf(:,sys_ind,flt_ind,t) - ...
%                 sum(Deltaf_Bu(:,:,sys_ind,flt_ind,t),2)- d(sys_ind,flt_ind,t)* fMode(flt_ind).g ==zeros(n,1)];
%              Constraint = [Constraint, Mode(sys_ind).C * z_c(:,sys_ind,flt_ind,t) + ...
%                 sum(Delta_Cx(:,:,sys_ind,flt_ind,t),2) + theta(:,sys_ind,flt_ind,t) == ...
%                 fMode(flt_ind).C * zf_c(:,sys_ind,flt_ind,t) + ...
%                 sum(Deltaf_Cx(:,:,sys_ind,flt_ind,t),2) + thetaf(:,sys_ind,flt_ind,t)];
%         end
%     end
% end


for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for flt_ind = 1: nf_mode  % fault mode index
            
             Constraint = [Constraint z_f(:,sys_ind,flt_ind,t) - Mode(sys_ind).A * z_c(:,sys_ind,flt_ind,t) - Mode(sys_ind).B * tt(:,sys_ind,flt_ind,t) -  d(sys_ind,flt_ind,t)* Mode(sys_ind).f == Delta(:,sys_ind)*Deltaf(sys_ind,flt_ind,t)];
             
             Constraint = [Constraint zf_f(:,sys_ind,flt_ind,t) - fMode(flt_ind).A * zf_c(:,sys_ind,flt_ind,t) - fMode(flt_ind).B * tt(:,sys_ind,flt_ind,t) - d(sys_ind,flt_ind,t)* (fMode(flt_ind).f) == Delta_f(:,flt_ind)*Deltaf_f(sys_ind,flt_ind,t) ];
            
             Constraint = [Constraint Mode(sys_ind).C * z_f(:,sys_ind,flt_ind,t) + theta(:,sys_ind,flt_ind,t)+ Mode(sys_ind).g == fMode(flt_ind).C * zf_f(:,sys_ind,flt_ind,t) + thetaf(:,sys_ind,flt_ind,t)+ fMode(flt_ind).g];
        end
    end
end

% for t = W+1:T
    Constraint = [Constraint sum(sum(sum(d(:,[3 4],1))))==1];
% end


% integer constraints

for t = 1:T
    Constraint = [Constraint sum(sum(d(:,:,t)))==1];
end

% current and future state constraints

for t = 1:T-1
    for flt_ind = 1:nf_mode
        Constraint = [Constraint sum(z_c(:,:,flt_ind,t+1),2)==sum(z_f(:,:,flt_ind,t),2)];
    end

    for sys_ind = 1:n_mode
        Constraint = [Constraint sum(zf_c(:,sys_ind,:,t+1),3)==sum(zf_f(:,sys_ind,:,t),3)];
    end
end
% admissible set constraints

for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for flt_ind = 1: nf_mode  % fault mode index
        for i = 1:n
            Constraint = [Constraint  d(sys_ind,flt_ind,t)* (M_l)<=z_c(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint z_c(i,sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)*M_u];
            Constraint = [Constraint  d(sys_ind,flt_ind,t)* (M_l)<=z_f(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint  z_f(i,sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)*M_u];
            Constraint = [Constraint  d(sys_ind,flt_ind,t)* (Mf_l)<=zf_c(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint  zf_c(i,sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)*Mf_u];
            Constraint = [Constraint  d(sys_ind,flt_ind,t)* (Mf_l)<=zf_f(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint  zf_f(i,sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)*Mf_u];
        end
        
        Constraint = [Constraint  d(sys_ind,flt_ind,t)* (-delta)<=Deltaf(sys_ind,flt_ind,t)];
        Constraint = [Constraint  Deltaf(sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)*delta];
        Constraint = [Constraint  d(sys_ind,flt_ind,t)* (-deltaf)<=Deltaf_f(sys_ind,flt_ind,t)];
        Constraint = [Constraint  Deltaf_f(sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)*deltaf];
        
        for i = 1:n_u
            Constraint = [Constraint  d(sys_ind,flt_ind,t)* (U_l)<=tt(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint tt(i,sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)*U_u];
        end
        
        for i = 1:n_y
            Constraint = [Constraint  d(sys_ind,flt_ind,t)* (-eps)<=theta(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint theta(i,sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)* eps];
            Constraint = [Constraint  d(sys_ind,flt_ind,t)* (-epsf)<=thetaf(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint thetaf(i,sys_ind,flt_ind,t)<= d(sys_ind,flt_ind,t)* epsf];
        end
        end 
    end
end


%% Setting up Options
options = sdpsettings('verbose',1,'solver',solver);

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