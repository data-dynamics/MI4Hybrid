function [Decision,sol] = Tdetect_uswa_milp_modified(SYS,SYS_f,T,N,Nf ...
    ,mn_bnd,pn_bnd,input_low, input_high, state_low, state_high,solver)

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
%   SYS_f -- a user-defined fault class (see StateSpace.m)
%   N -- user-defined normalization matrices for system model, which is 
%   itself a StateSpace class object(see StateSpace.m)
%   Nf -- user-defined normalization matrices for fault model, which is 
%   itself a StateSpace class object(see StateSpace.m)
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

%% If no normalization matrices are given for uncertainty


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
if(isempty(N))
    N_A = zeros(n,n,n_mode);
    N_B = zeros(n,n_u,n_mode);
    N_C = zeros(n_y,n,n_mode);
    N_D = zeros(n_y,n_u,n_mode);
    N_f = zeros(n,n_mode);
    N_g = zeros(n_y,n_mode);
    N = StateSpace(N_A,N_B,N_C,N_D,N_f,N_g);
end
if(isempty(Nf))
    Nf_A = zeros(n,n,nf_mode);
    Nf_B = zeros(n,n_u,nf_mode);
    Nf_C = zeros(n_y,n,nf_mode);
    Nf_D = zeros(n_y,n_u,nf_mode);
    Nf_f = zeros(n,nf_mode);
    Nf_g = zeros(n_y,nf_mode);
    Nf = StateSpace(Nf_A,Nf_B,Nf_C,Nf_D,Nf_f,Nf_g);
end
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

%% Check uncertainty normalization matrices
if(strcmp(N.mark,'ss')~=1&&strcmp(N.mark,'swss')~=1)
    error('The normalization matrices of system model must be a state-space model.');
end

if(strcmp(Nf.mark,'ss')~=1&&strcmp(Nf.mark,'swss')~=1)
    error('The normalization matrices of fault model must be a state-space model.');
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
        warning(['State upper bound is a scalar, converted to a vector with '...
        'identical entries.']);
    elseif (size(state_high,2)==2)
        state_high=[ones(n,1)*state_high(1,1) ones(n,1)*state_high(1,2)];
        warning(['State bound is a scalar, converted to a vector with '...
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
% Uncertainty normalization matrices
    Mode(i).N_A = N.mode(i).A;
    N_A(:,:,i) = N.mode(i).A;
    Mode(i).N_B = N.mode(i).B;
    N_B(:,:,i) = N.mode(i).B;
    Mode(i).N_C = N.mode(i).C;
    N_C(:,:,i) = N.mode(i).C;
    Mode(i).N_D = N.mode(i).D;
    N_D(:,:,i) = N.mode(i).D;
    Mode(i).N_f = N.mode(i).f;
    N_f(:,i) = N.mode(i).f;
    Mode(i).N_g = N.mode(i).g;
    N_g(:,i) = N.mode(i).g;
end

%% Initiate fault modes
for i = 1: nf_mode
    fMode(i).A = SYS_f.mode(i).A;
    fMode(i).B = SYS_f.mode(i).B;
    fMode(i).C = SYS_f.mode(i).C;
    fMode(i).D = SYS_f.mode(i).D;
    fMode(i).f = SYS_f.mode(i).f;
    fMode(i).g = SYS_f.mode(i).g;
% Uncertainty normalization matrices
    fMode(i).N_A = Nf.mode(i).A;
    Nf_A(:,:,i) = Nf.mode(i).A;
    fMode(i).N_B = Nf.mode(i).B;
    Nf_B(:,:,i) = Nf.mode(i).B;
    fMode(i).N_C = Nf.mode(i).C;
    Nf_C(:,:,i) = Nf.mode(i).C;
    fMode(i).N_D = Nf.mode(i).D;
    Nf_D(:,:,i) = Nf.mode(i).D;
    fMode(i).N_f = Nf.mode(i).f;
    Nf_f(:,i) = Nf.mode(i).f;
    fMode(i).N_g = Nf.mode(i).g;
    Nf_g(:,i) = Nf.mode(i).g;
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
% binary variables
d = binvar(n_mode,nf_mode,T,'full');      
d_A = binvar(n,n,4,n_mode,nf_mode,T,'full');
d_B = binvar(n,n_u,4,n_mode,nf_mode,T,'full');
d_C = binvar(n_y,n,4,n_mode,nf_mode,T,'full');
d_D = binvar(n_y,n_u,4,n_mode,nf_mode,T,'full');
fd_A = binvar(n,n,4,n_mode,nf_mode,T,'full');
fd_B = binvar(n,n_u,4,n_mode,nf_mode,T,'full');
fd_C = binvar(n_y,n,4,n_mode,nf_mode,T,'full');
fd_D = binvar(n_y,n_u,4,n_mode,nf_mode,T,'full');
% real variables
cz = sdpvar(n,n_mode,nf_mode,T,'full');  % current state
fz = sdpvar(n,n_mode,nf_mode,T,'full');  % future state
cz_f = sdpvar(n,n_mode,nf_mode,T,'full');  % fault current state
fz_f = sdpvar(n,n_mode,nf_mode,T,'full');  % fault future state
theta = sdpvar(n_y,n_mode,nf_mode,T,'full'); % measurement noise (system)
thetaf = sdpvar(n_y,n_mode,nf_mode,T,'full'); % measurement noise (fault)
eta = sdpvar(n,n_mode,nf_mode,T,'full'); % process noise (system)
etaf = sdpvar(n,n_mode,nf_mode,T,'full'); % process noise (fault)
tt = sdpvar(n_u,n_mode,nf_mode,T,'full'); % input 

% Uncertainty variables for system model
DA = sdpvar(n,n,n_mode,nf_mode,T,'full');
DB = sdpvar(n,n_u,n_mode,nf_mode,T,'full');
DC = sdpvar(n_y,n,n_mode,nf_mode,T,'full');
DD = sdpvar(n_y,n_u,n_mode,nf_mode,T,'full');
Df = sdpvar(n,n_mode,nf_mode,T,'full');
Dg = sdpvar(n_y,n_mode,nf_mode,T,'full');
SA = sdpvar(n,n,4,n_mode,nf_mode,T,'full');
SB = sdpvar(n,n_u,4,n_mode,nf_mode,T,'full');
SC = sdpvar(n_y,n,4,n_mode,nf_mode,T,'full');
SD = sdpvar(n_y,n_u,4,n_mode,nf_mode,T,'full');
zc_A = sdpvar(n,4,n_mode,nf_mode,T,'full');
zc_C = sdpvar(n,4,n_mode,nf_mode,T,'full');
tt_B = sdpvar(n_u,4,n_mode,nf_mode,T,'full');
tt_D = sdpvar(n_u,4,n_mode,nf_mode,T,'full');

% Uncertainty variables for fault model
fDA = sdpvar(n,n,n_mode,nf_mode,T,'full');
fDB = sdpvar(n,n_u,n_mode,nf_mode,T,'full');
fDC = sdpvar(n_y,n,n_mode,nf_mode,T,'full');
fDD = sdpvar(n_y,n_u,n_mode,nf_mode,T,'full');
fDf = sdpvar(n,n_mode,nf_mode,T,'full');
fDg = sdpvar(n_y,n_mode,nf_mode,T,'full');
fSA = sdpvar(n,n,4,n_mode,nf_mode,T,'full');
fSB = sdpvar(n,n_u,4,n_mode,nf_mode,T,'full');
fSC = sdpvar(n_y,n,4,n_mode,nf_mode,T,'full');
fSD = sdpvar(n_y,n_u,4,n_mode,nf_mode,T,'full');
fzc_A = sdpvar(n,4,n_mode,nf_mode,T,'full');
fzc_C = sdpvar(n,4,n_mode,nf_mode,T,'full');
ftt_B = sdpvar(n_u,4,n_mode,nf_mode,T,'full');
ftt_D = sdpvar(n_u,4,n_mode,nf_mode,T,'full');

%% Creating Constraints
Constraint = [];

% state and output equation constraints:
for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for flt_ind = 1: nf_mode  % fault mode index
            
         Constraint = [Constraint fz(:,sys_ind,flt_ind,t) - ...
             Mode(sys_ind).A * cz(:,sys_ind,flt_ind,t) - ...
             sum(DA(:,:,sys_ind,flt_ind,t),2) - ...
             Mode(sys_ind).B * tt(:,sys_ind,flt_ind,t) -  ...
             sum(DB(:,:,sys_ind,flt_ind,t),2)-...
             eta(:,sys_ind,flt_ind,t)-Df(:,sys_ind,flt_ind,t) == ...
             d(sys_ind,flt_ind,t)* Mode(sys_ind).f];
             
         Constraint = [Constraint fz_f(:,sys_ind,flt_ind,t) - ...
             fMode(flt_ind).A * cz_f(:,sys_ind,flt_ind,t) - ...
             sum(fDA(:,:,sys_ind,flt_ind,t),2)...
             - fMode(flt_ind).B* tt(:,sys_ind,flt_ind,t)-  ...
             sum(fDB(:,:,sys_ind,flt_ind,t),2)-...
             etaf(:,sys_ind,flt_ind,t)-fDf(:,sys_ind,flt_ind,t) == ...
             d(sys_ind,flt_ind,t)* (fMode(flt_ind).f)];
            
         Constraint = [Constraint Mode(sys_ind).C * fz(:,sys_ind,flt_ind,t)...
             + sum(DC(:,:,sys_ind,flt_ind,t),2)+...
             Mode(sys_ind).D*tt(:,sys_ind,flt_ind,t)+...
             sum(DD(:,:,sys_ind,flt_ind,t),2)+...
             theta(:,sys_ind,flt_ind,t)+ d(sys_ind,flt_ind,t)*Mode(sys_ind).g+ ...
             Dg(:,sys_ind,flt_ind,t) == fMode(flt_ind).C *...
             fz_f(:,sys_ind,flt_ind,t) + sum(fDC(:,:,sys_ind,flt_ind,t),2)+...
             fMode(flt_ind).D*tt(:,sys_ind,flt_ind,t) + sum(fDD(:,:,sys_ind,flt_ind,t),2)+...
             thetaf(:,sys_ind,flt_ind,t)+ d(sys_ind,flt_ind,t)*fMode(flt_ind).g +...
             fDg(:,sys_ind,flt_ind,t)];
        end
    end
end

% integer constraints
for t = 1:T
    Constraint = [Constraint sum(sum(d(:,:,t)))==1];
end

% same initial condition constraint
Constraint = [Constraint sum(sum(cz(:,:,:,1),3),2)==...
    sum(sum(cz_f(:,:,:,1),3),2)];

% current and future state constraints
for t = 1:T-1
    for flt_ind = 1:nf_mode
        Constraint = [Constraint sum(cz(:,:,flt_ind,t+1),2)==...
            sum(fz(:,:,flt_ind,t),2)];
    end

    for sys_ind = 1:n_mode
        Constraint = [Constraint sum(cz_f(:,sys_ind,:,t+1),3)==...
            sum(fz_f(:,sys_ind,:,t),3)];
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
                cz(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint cz(i,sys_ind,flt_ind,t)<= ...
                d(sys_ind,flt_ind,t)*M_u(i)];
            Constraint = [Constraint d(sys_ind,flt_ind,t)* (Mf_l(i))<=...
                cz_f(i,sys_ind,flt_ind,t)];
            Constraint = [Constraint cz_f(i,sys_ind,flt_ind,t)<=...
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

%% Uncertainty constraints
%% f and g matrices
if (isequal(N_f,zeros(n,n_mode)))
    Constraint = [Constraint,Df == zeros(n,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for i = 1:n
                    Constraint = [Constraint  -d(sys_ind,flt_ind,t)*abs(Mode(sys_ind).N_f(i))<=...
                        Df(i,sys_ind,flt_ind,t)];
                    Constraint = [Constraint Df(i,sys_ind,flt_ind,t)<= ...
                        d(sys_ind,flt_ind,t)*abs(Mode(sys_ind).N_f(i))];
                end
                % Extra constraints, because all the uncertainty is caused by
                % Ambient temperature
                for i = 3:n-1
                    for j = i+1:n
                        Constraint = [Constraint, Df(i,sys_ind,flt_ind,t)...
                            * Mode(sys_ind).N_f(j)==Df(j,sys_ind,flt_ind,t)...
                            *Mode(sys_ind).N_f(i)];
                    end
                end
            end
        end
    end
end
if (isequal(Nf_f,zeros(n,nf_mode)));
    Constraint = [Constraint,fDf == zeros(n,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for i = 1:n
                    Constraint = [Constraint  -d(sys_ind,flt_ind,t)*abs(fMode(flt_ind).N_f(i))<=...
                        fDf(i,sys_ind,flt_ind,t)];
                    Constraint = [Constraint fDf(i,sys_ind,flt_ind,t)<= ...
                        d(sys_ind,flt_ind,t)*abs(fMode(flt_ind).N_f(i))];
                end
                % Extra constraints, because all the uncertainty is caused by
                % Ambient temperature
                for i = 3:n-1
                    for j = i+1:n
                        Constraint = [Constraint, fDf(i,sys_ind,flt_ind,t)...
                            * fMode(flt_ind).N_f(j)==fDf(j,sys_ind,flt_ind,t)...
                            fMode(flt_ind).N_f(i)];
                    end
                end
            end
        end
    end
end

if (isequal(N_g, zeros(n_y,n_mode)));
    Constraint = [Constraint,Dg == zeros(n_y,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for i = 1:n_y
                    Constraint = [Constraint  -d(sys_ind,flt_ind,t)*abs(Mode(sys_ind).N_g(i))<=...
                        Dg(i,sys_ind,flt_ind,t)];
                    Constraint = [Constraint Dg(i,sys_ind,flt_ind,t)<= ...
                        d(sys_ind,flt_ind,t)*abs(Mode(sys_ind).N_g(i))];
                end
            end
        end
    end
end

if (isequal(Nf_g, zeros(n_y,nf_mode)));
    Constraint = [Constraint,fDg == zeros(n_y,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for i = 1:n_y
                    Constraint = [Constraint  -d(sys_ind,flt_ind,t)*abs(fMode(flt_ind).N_g(i))<=...
                        fDg(i,sys_ind,flt_ind,t)];
                    Constraint = [Constraint fDg(i,sys_ind,flt_ind,t)<= ...
                        d(sys_ind,flt_ind,t)*abs(fMode(flt_ind).N_g(i))];
                end
            end
        end
    end
end


%% A,B,C,D matrices
%% System Uncertainty
% Uncertainty in A
if (isequal(N_A, zeros(n,n,n_mode)));
    Constraint = [Constraint,DA == zeros(n,n,n_mode,nf_mode,T)];
else
for t = 1:T    % time index
    for sys_ind = 1: n_mode  % system mode index
        for flt_ind = 1: nf_mode  % fault mode index
            for j = 1:n
                for i= 1:n
                    Constraint = [Constraint, SA(i,j,1,sys_ind,flt_ind,t)>=0, ...
                        zc_A(j,1,sys_ind,flt_ind,t)>=0, ...
                        SA(i,j,1,sys_ind,flt_ind,t)-abs(Mode(sys_ind).N_A(i,j))*...
                        zc_A(j,1,sys_ind,flt_ind,t)<=0];
                    Constraint = [Constraint, SA(i,j,2,sys_ind,flt_ind,t)<=0, ...
                        zc_A(j,2,sys_ind,flt_ind,t)<=0, ...
                        SA(i,j,2,sys_ind,flt_ind,t)-abs(Mode(sys_ind).N_A(i,j))*...
                        zc_A(j,2,sys_ind,flt_ind,t)>=0];
                    Constraint = [Constraint, SA(i,j,3,sys_ind,flt_ind,t)>=0, ...
                        zc_A(j,3,sys_ind,flt_ind,t)<=0, ...
                        SA(i,j,3,sys_ind,flt_ind,t)+abs(Mode(sys_ind).N_A(i,j))*...
                        zc_A(j,3,sys_ind,flt_ind,t)<=0];
                    Constraint = [Constraint, SA(i,j,4,sys_ind,flt_ind,t)<=0, ...
                        zc_A(j,4,sys_ind,flt_ind,t)>=0, ...
                        SA(i,j,4,sys_ind,flt_ind,t)+abs(Mode(sys_ind).N_A(i,j))*...
                        zc_A(j,4,sys_ind,flt_ind,t)>=0];
                    for q=1:4
                        Constraint = [Constraint, SA(i,j,q,sys_ind,flt_ind,t)<=...
                            d_A(i,j,q,sys_ind,flt_ind,t)*M_u(j)*abs(Mode(sys_ind).N_A(i,j)), ...
                            SA(i,j,q,sys_ind,flt_ind,t)>=...
                            d_A(i,j,q,sys_ind,flt_ind,t)*M_l(j)*abs(Mode(sys_ind).N_A(i,j))];
                    end
                    Constraint = [Constraint, sum(SA(i,j,:,sys_ind,flt_ind,t)) == DA(i,j,sys_ind,flt_ind,t)];
                    Constraint = [Constraint, sum(d_A(i,j,:,sys_ind,flt_ind,t)) == d(sys_ind,flt_ind,t)];
                end
                Constraint = [Constraint, sum(zc_A(j,:,sys_ind,flt_ind,t))>= d(sys_ind,flt_ind,t)*M_l(j)];
                Constraint = [Constraint, sum(zc_A(j,:,sys_ind,flt_ind,t))<= d(sys_ind,flt_ind,t)*M_u(j)];
                Constraint = [Constraint, sum(zc_A(j,:,sys_ind,flt_ind,t)) == cz(j,sys_ind,flt_ind,t)]; 
            end
        end
    end
end
end
% Uncertainty in B
if (isequal(N_B, zeros(n,n_u,n_mode)));
    Constraint = [Constraint, DB == zeros(n,n_u,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                
                for j = 1:n_u
                    for i= 1:n
                        Constraint = [Constraint, SB(i,j,1,sys_ind,flt_ind,t)>=0, ...
                            tt_B(j,1,sys_ind,flt_ind,t)>=0, ...
                            SB(i,j,1,sys_ind,flt_ind,t)-abs(Mode(sys_ind).N_B(i,j))*...
                            tt_B(j,1,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, SB(i,j,2,sys_ind,flt_ind,t)<=0, ...
                            tt_B(j,2,sys_ind,flt_ind,t)<=0, ...
                            SB(i,j,2,sys_ind,flt_ind,t)-abs(Mode(sys_ind).N_B(i,j))*...
                            tt_B(j,2,sys_ind,flt_ind,t)>=0];
                        Constraint = [Constraint, SB(i,j,3,sys_ind,flt_ind,t)>=0, ...
                            tt_B(j,3,sys_ind,flt_ind,t)<=0, ...
                            SB(i,j,3,sys_ind,flt_ind,t)+abs(Mode(sys_ind).N_B(i,j))*...
                            tt_B(j,3,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, SB(i,j,4,sys_ind,flt_ind,t)<=0, ...
                            tt_B(j,4,sys_ind,flt_ind,t)>=0, ...
                            SB(i,j,4,sys_ind,flt_ind,t)+abs(Mode(sys_ind).N_B(i,j))*...
                            tt_B(j,4,sys_ind,flt_ind,t)>=0];
                        for q=1:4
                            Constraint = [Constraint, SB(i,j,q,sys_ind,flt_ind,t)<=...
                                d_B(i,j,q,sys_ind,flt_ind,t)*abs(Mode(sys_ind).N_B(i,j))*U_u(j), ...
                                SB(i,j,q,sys_ind,flt_ind,t)>=...
                                d_B(i,j,q,sys_ind,flt_ind,t)*abs(Mode(sys_ind).N_B(i,j))*U_l(j)];
                        end
                        Constraint = [Constraint, sum(SB(i,j,:,sys_ind,flt_ind,t)) == DB(i,j,sys_ind,flt_ind,t)];
                        Constraint = [Constraint, sum(d_B(i,j,:,sys_ind,flt_ind,t)) == d(sys_ind,flt_ind,t)];
                    end
                    Constraint = [Constraint, sum(tt_B(j,:,sys_ind,flt_ind,t))>= d(sys_ind,flt_ind,t)*U_l(j)];
                    Constraint = [Constraint, sum(tt_B(j,:,sys_ind,flt_ind,t))<= d(sys_ind,flt_ind,t)*U_u(j)];
                    Constraint = [Constraint, sum(tt_B(j,:,sys_ind,flt_ind,t)) == tt(j,sys_ind,flt_ind,t)];
                end
            end
        end
    end
end
% Uncertainty for C
if (isequal(N_C, zeros(n_y,n,n_mode)));
    Constraint = [Constraint, DC == zeros(n_y,n,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for j = 1:n
                    for i= 1:n_y
                        Constraint = [Constraint, SC(i,j,1,sys_ind,flt_ind,t)>=0, ...
                            zc_C(j,1,sys_ind,flt_ind,t)>=0, ...
                            SC(i,j,1,sys_ind,flt_ind,t)-abs(Mode(sys_ind).N_C(i,j))*...
                            zc_C(j,1,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, SC(i,j,2,sys_ind,flt_ind,t)<=0, ...
                            zc_C(j,2,sys_ind,flt_ind,t)<=0, ...
                            SC(i,j,2,sys_ind,flt_ind,t)-abs(Mode(sys_ind).N_C(i,j))*...
                            zc_C(j,2,sys_ind,flt_ind,t)>=0];
                        Constraint = [Constraint, SC(i,j,3,sys_ind,flt_ind,t)>=0, ...
                            zc_C(j,3,sys_ind,flt_ind,t)<=0, ...
                            SC(i,j,3,sys_ind,flt_ind,t)+abs(Mode(sys_ind).N_C(i,j))*...
                            zc_C(j,3,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, SC(i,j,4,sys_ind,flt_ind,t)<=0, ...
                            zc_C(j,4,sys_ind,flt_ind,t)>=0, ...
                            SC(i,j,4,sys_ind,flt_ind,t)+abs(Mode(sys_ind).N_C(i,j))*...
                            zc_C(j,4,sys_ind,flt_ind,t)>=0];
                        for q=1:4
                            Constraint = [Constraint, SC(i,j,q,sys_ind,flt_ind,t)<=...
                                d_C(i,j,q,sys_ind,flt_ind,t)*M_u(j)*abs(Mode(sys_ind).N_C(i,j)), ...
                                SC(i,j,q,sys_ind,flt_ind,t)>=...
                                d_C(i,j,q,sys_ind,flt_ind,t)*M_l(j)*abs(Mode(sys_ind).N_C(i,j))];
                        end
                        Constraint = [Constraint, sum(d_C(i,j,:,sys_ind,flt_ind,t)) == d(sys_ind,flt_ind,t)];
                        Constraint = [Constraint, sum(SC(i,j,:,sys_ind,flt_ind,t)) == DC(i,j,sys_ind,flt_ind,t)];
                        
                    end
                    Constraint = [Constraint, sum(zc_C(j,:,sys_ind,flt_ind,t))>= d(sys_ind,flt_ind,t)*M_l(j)];
                    Constraint = [Constraint, sum(zc_C(j,:,sys_ind,flt_ind,t))<= d(sys_ind,flt_ind,t)*M_u(j)];
                    Constraint = [Constraint, sum(zc_C(j,:,sys_ind,flt_ind,t)) == fz(j,sys_ind,flt_ind,t)];
                end
            end
        end
    end
end

% Uncertinty for D
if (isequal(N_D, zeros(n_y,n_u,n_mode)));
    Constraint = [Constraint, DD == zeros(n_y,n_u,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for j = 1:n_u
                    for i= 1:n_y
                        Constraint = [Constraint, SD(i,j,1,sys_ind,flt_ind,t)>=0, ...
                            tt_D(j,1,sys_ind,flt_ind,t)>=0, ...
                            SD(i,j,1,sys_ind,flt_ind,t)-abs(Mode(sys_ind).N_D(i,j))*...
                            tt_D(j,1,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, SD(i,j,2,sys_ind,flt_ind,t)<=0, ...
                            tt_D(j,2,sys_ind,flt_ind,t)<=0, ...
                            SD(i,j,2,sys_ind,flt_ind,t)-abs(Mode(sys_ind).N_D(i,j))*...
                            tt_D(j,2,sys_ind,flt_ind,t)>=0];
                        Constraint = [Constraint, SD(i,j,3,sys_ind,flt_ind,t)>=0, ...
                            tt_D(j,3,sys_ind,flt_ind,t)<=0, ...
                            SD(i,j,3,sys_ind,flt_ind,t)+abs(Mode(sys_ind).N_D(i,j))*...
                            tt_D(j,3,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, SD(i,j,4,sys_ind,flt_ind,t)<=0, ...
                            tt_D(j,4,sys_ind,flt_ind,t)>=0, ...
                            SD(i,j,4,sys_ind,flt_ind,t)+abs(Mode(sys_ind).N_D(i,j))*...
                            tt_D(j,4,sys_ind,flt_ind,t)>=0];
                        for q=1:4
                            Constraint = [Constraint, SD(i,j,q,sys_ind,flt_ind,t)<=...
                                d_D(i,j,q,sys_ind,flt_ind,t)*abs(Mode(sys_ind).N_D(i,j))*U_u(j), ...
                                SD(i,j,q,sys_ind,flt_ind,t)>=...
                                d_D(i,j,q,sys_ind,flt_ind,t)*abs(Mode(sys_ind).N_D(i,j))*U_l(j)];
                        end
                        Constraint = [Constraint, sum(SD(i,j,:,sys_ind,flt_ind,t)) == DD(i,j,sys_ind,flt_ind,t)];
                        Constraint = [Constraint, sum(d_D(i,j,:,sys_ind,flt_ind,t)) == d(sys_ind,flt_ind,t)];
                    end
                    Constraint = [Constraint, sum(tt_D(j,:,sys_ind,flt_ind,t))>= d(sys_ind,flt_ind,t)*U_l(j)];
                    Constraint = [Constraint, sum(tt_D(j,:,sys_ind,flt_ind,t))<= d(sys_ind,flt_ind,t)*U_u(j)];
                    Constraint = [Constraint, sum(tt_D(j,:,sys_ind,flt_ind,t)) == tt(j,sys_ind,flt_ind,t)];
                end
            end
        end
    end
end

%% Fault Uncertainty
% Uncertainty in A
if (isequal(Nf_A, zeros(n,n,nf_mode)));
    Constraint = [Constraint,fDA == zeros(n,n,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for j = 1:n
                    for i= 1:n
                        Constraint = [Constraint, fSA(i,j,1,sys_ind,flt_ind,t)>=0, ...
                            fzc_A(j,1,sys_ind,flt_ind,t)>=0, ...
                            fSA(i,j,1,sys_ind,flt_ind,t)-abs(fMode(flt_ind).N_A(i,j))*...
                            fzc_A(j,1,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, fSA(i,j,2,sys_ind,flt_ind,t)<=0, ...
                            fzc_A(j,2,sys_ind,flt_ind,t)<=0, ...
                            fSA(i,j,2,sys_ind,flt_ind,t)-abs(fMode(flt_ind).N_A(i,j))*...
                            fzc_A(j,2,sys_ind,flt_ind,t)>=0];
                        Constraint = [Constraint, fSA(i,j,3,sys_ind,flt_ind,t)>=0, ...
                            fzc_A(j,3,sys_ind,flt_ind,t)<=0, ...
                            fSA(i,j,3,sys_ind,flt_ind,t)+abs(fMode(flt_ind).N_A(i,j))*...
                            fzc_A(j,3,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, fSA(i,j,4,sys_ind,flt_ind,t)<=0, ...
                            fzc_A(j,4,sys_ind,flt_ind,t)>=0, ...
                            fSA(i,j,4,sys_ind,flt_ind,t)+abs(fMode(flt_ind).N_A(i,j))*...
                            fzc_A(j,4,sys_ind,flt_ind,t)>=0];
                        for q=1:4
                            Constraint = [Constraint, fSA(i,j,q,sys_ind,flt_ind,t)<=...
                                fd_A(i,j,q,sys_ind,flt_ind,t)*Mf_u(j)*abs(fMode(flt_ind).N_A(i,j)), ...
                                fSA(i,j,q,sys_ind,flt_ind,t)>=...
                                fd_A(i,j,q,sys_ind,flt_ind,t)*Mf_l(j)*abs(fMode(flt_ind).N_A(i,j))];
                        end
                        Constraint = [Constraint, sum(fSA(i,j,:,sys_ind,flt_ind,t)) == fDA(i,j,sys_ind,flt_ind,t)];
                        Constraint = [Constraint, sum(fd_A(i,j,:,sys_ind,flt_ind,t)) == d(sys_ind,flt_ind,t)];
                    end
                    Constraint = [Constraint, sum(fzc_A(j,:,sys_ind,flt_ind,t))>= d(sys_ind,flt_ind,t)*Mf_l(j)];
                    Constraint = [Constraint, sum(fzc_A(j,:,sys_ind,flt_ind,t))<= d(sys_ind,flt_ind,t)*Mf_u(j)];
                    Constraint = [Constraint, sum(fzc_A(j,:,sys_ind,flt_ind,t)) == cz_f(j,sys_ind,flt_ind,t)];
                end
            end
        end
    end
end
% Uncertainty in B
if (isequal(Nf_B, zeros(n,n_u,nf_mode)));
    Constraint = [Constraint,fDB == zeros(n,n_u,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for j = 1:n_u
                    for i= 1:n
                        Constraint = [Constraint, fSB(i,j,1,sys_ind,flt_ind,t)>=0, ...
                            ftt_B(j,1,sys_ind,flt_ind,t)>=0, ...
                            fSB(i,j,1,sys_ind,flt_ind,t)-abs(fMode(flt_ind).N_B(i,j))*...
                            ftt_B(j,1,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, fSB(i,j,2,sys_ind,flt_ind,t)<=0, ...
                            ftt_B(j,2,sys_ind,flt_ind,t)<=0, ...
                            fSB(i,j,2,sys_ind,flt_ind,t)-abs(fMode(flt_ind).N_B(i,j))*...
                            ftt_B(j,2,sys_ind,flt_ind,t)>=0];
                        Constraint = [Constraint, fSB(i,j,3,sys_ind,flt_ind,t)>=0, ...
                            ftt_B(j,3,sys_ind,flt_ind,t)<=0, ...
                            fSB(i,j,3,sys_ind,flt_ind,t)+abs(fMode(flt_ind).N_B(i,j))*...
                            ftt_B(j,3,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, fSB(i,j,4,sys_ind,flt_ind,t)<=0, ...
                            ftt_B(j,4,sys_ind,flt_ind,t)>=0, ...
                            fSB(i,j,4,sys_ind,flt_ind,t)+abs(fMode(flt_ind).N_B(i,j))*...
                            ftt_B(j,4,sys_ind,flt_ind,t)>=0];
                        for q=1:4
                            Constraint = [Constraint, fSB(i,j,q,sys_ind,flt_ind,t)<=...
                                d_B(i,j,q,sys_ind,flt_ind,t)*abs(fMode(flt_ind).N_B(i,j))*U_u(j), ...
                                SB(i,j,q,sys_ind,flt_ind,t)>=...
                                d_B(i,j,q,sys_ind,flt_ind,t)*abs(fMode(flt_ind).N_B(i,j))*U_l(j)];
                        end
                        Constraint = [Constraint, sum(SB(i,j,:,sys_ind,flt_ind,t)) == DB(i,j,sys_ind,flt_ind,t)];
                        Constraint = [Constraint, sum(d_B(i,j,:,sys_ind,flt_ind,t)) == d(sys_ind,flt_ind,t)];
                    end
                    Constraint = [Constraint, sum(tt_B(j,:,sys_ind,flt_ind,t))>= d(sys_ind,flt_ind,t)*U_l(j)];
                    Constraint = [Constraint, sum(tt_B(j,:,sys_ind,flt_ind,t))<= d(sys_ind,flt_ind,t)*U_u(j)];
                    Constraint = [Constraint, sum(tt_B(j,:,sys_ind,flt_ind,t)) == tt(j,sys_ind,flt_ind,t)];
                end
            end
        end
    end
end
% Uncertainty for C
if (isequal(Nf_C, zeros(n_y,n,nf_mode)));
    Constraint = [Constraint,fDC == zeros(n_y,n,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for j = 1:n
                    for i= 1:n_y
                        Constraint = [Constraint, fSC(i,j,1,sys_ind,flt_ind,t)>=0, ...
                            fzc_C(j,1,sys_ind,flt_ind,t)>=0, ...
                            fSC(i,j,1,sys_ind,flt_ind,t)-abs(fMode(flt_ind).N_C(i,j))*...
                            fzc_C(j,1,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, fSC(i,j,2,sys_ind,flt_ind,t)<=0, ...
                            fzc_C(j,2,sys_ind,flt_ind,t)<=0, ...
                            fSC(i,j,2,sys_ind,flt_ind,t)-abs(fMode(flt_ind).N_C(i,j))*...
                            fzc_C(j,2,sys_ind,flt_ind,t)>=0];
                        Constraint = [Constraint, fSC(i,j,3,sys_ind,flt_ind,t)>=0, ...
                            fzc_C(j,3,sys_ind,flt_ind,t)<=0, ...
                            fSC(i,j,3,sys_ind,flt_ind,t)+abs(fMode(flt_ind).N_C(i,j))*...
                            fzc_C(j,3,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, SC(i,j,4,sys_ind,flt_ind,t)<=0, ...
                            fzc_C(j,4,sys_ind,flt_ind,t)>=0, ...
                            fSC(i,j,4,sys_ind,flt_ind,t)+abs(fMode(flt_ind).N_C(i,j))*...
                            fzc_C(j,4,sys_ind,flt_ind,t)>=0];
                        for q=1:4
                            Constraint = [Constraint, fSC(i,j,q,sys_ind,flt_ind,t)<=...
                                fd_C(i,j,q,sys_ind,flt_ind,t)*Mf_u(j)*abs(fMode(flt_ind).N_C(i,j)), ...
                                fSC(i,j,q,sys_ind,flt_ind,t)>=...
                                fd_C(i,j,q,sys_ind,flt_ind,t)*Mf_l(j)*abs(fMode(flt_ind).N_C(i,j))];
                        end
                        Constraint = [Constraint, sum(fd_C(i,j,:,sys_ind,flt_ind,t)) == d(sys_ind,flt_ind,t)];
                        Constraint = [Constraint, sum(fSC(i,j,:,sys_ind,flt_ind,t)) == fDC(i,j,sys_ind,flt_ind,t)];
                        
                    end
                    Constraint = [Constraint, sum(fzc_C(j,:,sys_ind,flt_ind,t))>= d(sys_ind,flt_ind,t)*Mf_l(j)];
                    Constraint = [Constraint, sum(fzc_C(j,:,sys_ind,flt_ind,t))<= d(sys_ind,flt_ind,t)*Mf_u(j)];
                    Constraint = [Constraint, sum(fzc_C(j,:,sys_ind,flt_ind,t)) == fz_f(j,sys_ind,flt_ind,t)];
                end
            end
        end
    end
end
% Uncertinty for D
if (isequal(Nf_D, zeros(n_y,n_u,nf_mode)));
    Constraint = [Constraint,fDD == zeros(n_y,n_u,n_mode,nf_mode,T)];
else
    for t = 1:T    % time index
        for sys_ind = 1: n_mode  % system mode index
            for flt_ind = 1: nf_mode  % fault mode index
                for j = 1:n_u
                    for i= 1:n_y
                        Constraint = [Constraint, fSD(i,j,1,sys_ind,flt_ind,t)>=0, ...
                            ftt_D(j,1,sys_ind,flt_ind,t)>=0, ...
                            fSD(i,j,1,sys_ind,flt_ind,t)-abs(fMode(flt_ind).N_D(i,j))*...
                            ftt_D(j,1,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, fSD(i,j,2,sys_ind,flt_ind,t)<=0, ...
                            ftt_D(j,2,sys_ind,flt_ind,t)<=0, ...
                            fSD(i,j,2,sys_ind,flt_ind,t)-abs(fMode(flt_ind).N_D(i,j))*...
                            ftt_D(j,2,sys_ind,flt_ind,t)>=0];
                        Constraint = [Constraint, fSD(i,j,3,sys_ind,flt_ind,t)>=0, ...
                            ftt_D(j,3,sys_ind,flt_ind,t)<=0, ...
                            fSD(i,j,3,sys_ind,flt_ind,t)+abs(fMode(flt_ind).N_D(i,j))*...
                            ftt_D(j,3,sys_ind,flt_ind,t)<=0];
                        Constraint = [Constraint, fSD(i,j,4,sys_ind,flt_ind,t)<=0, ...
                            ftt_D(j,4,sys_ind,flt_ind,t)>=0, ...
                            fSD(i,j,4,sys_ind,flt_ind,t)+abs(fMode(flt_ind).N_D(i,j))*...
                            ftt_D(j,4,sys_ind,flt_ind,t)>=0];
                        for q=1:4
                            Constraint = [Constraint, fSD(i,j,q,sys_ind,flt_ind,t)<=...
                                fd_D(i,j,q,sys_ind,flt_ind,t)*abs(fMode(flt_ind).N_D(i,j))*U_u(j), ...
                                fSD(i,j,q,sys_ind,flt_ind,t)>=...
                                fd_D(i,j,q,sys_ind,flt_ind,t)*abs(fMode(flt_ind).N_D(i,j))*U_l(j)];
                        end
                        Constraint = [Constraint, sum(fSD(i,j,:,sys_ind,flt_ind,t)) == fDD(i,j,sys_ind,flt_ind,t)];
                        Constraint = [Constraint, sum(fd_D(i,j,:,sys_ind,flt_ind,t)) == d(sys_ind,flt_ind,t)];
                    end
                    Constraint = [Constraint, sum(ftt_D(j,:,sys_ind,flt_ind,t))>= d(sys_ind,flt_ind,t)*U_l(j)];
                    Constraint = [Constraint, sum(ftt_D(j,:,sys_ind,flt_ind,t))<= d(sys_ind,flt_ind,t)*U_u(j)];
                    Constraint = [Constraint, sum(ftt_D(j,:,sys_ind,flt_ind,t)) == tt(j,sys_ind,flt_ind,t)];
                end
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