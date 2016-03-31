function [miflag, EPS]=invalidation_poly(sys,u,y,BP, BM, LS, US, param, ...
    verbose)

% This function uses input/output data and user specified noise/state bound
% (and uncertainty bound, if any) to determine if a non-switched polynomial
% model is invalidated. The returned value can also indicate how "far" the
% system is from being not invalidated.
%
% Input:
%   sys -- a user defined polynomial model (see PolyModel.m or UnPolyModel.m)
%   u -- an ni-by-T input sequence to the polynomial model where ni is the
%        input dimension and T is the time horizon
%   y -- an n-by-T output sequence to the polynomial model where n is the
%        output/state dimension and T is the time horizon
%   BP -- an n-D vector specifying the bounds on process noise
%   BM -- an n-D vector specifying the bounds on measurement noise
%   LS -- an n-D vector specifying the lower bounds on the states
%   US -- an n-D vector specifying the upper bounds on the states
%   param -- User can refine the solution obtained from the SDP relaxation
%            by setting the parameter:
%               param.POPsolver = 'active-set';
%               param.POPsolver = 'interior-point';
%               param.POPsolver = 'trust-region-reflective';
%               param.POPsolver = 'lsqnonlin';
%            See SparsePOP userGuide.pdf for more information about
%               param.relaxOrder, param.eqTolerance, and param.scalingSW
%   verbose -- an integer (varying from 0 to 2) used to control how much
%              solving information to be printed
%
% Output:
%   miflag -- a flag indicating whether the model is invalidated (0 for not
%             invalidated and 1 for invalidated)
%   EPS -- the optimal value for the objective function. If this value is
%            NaN, it means that the optimization problem is infeasible (and
%            the model is invalidated).
%
% Syntax:
%   [miflag,sigma]=invalidation_poly(sys,u,y,BP,BM);
%   [miflag,sigma]=invalidation_poly(sys,u,y,BP,BM,LS,US);
%   [miflag,sigma]=invalidation_poly(sys,u,y,BP,BM,LS,US,param);
%   [miflag,sigma]=invalidation_poly(sys,u,y,BP,BM,LS,US,param,verbose);
%
% Note: the objective function can be different for different cases (e.g.
% with or without uncertainty). Details are discussed in the paper "Model
% (In)validation and Fault Detection for Systems with Polynomial
% State-Space Models" by F. Harirchi, Z. Luo, and N. Ozay. However, here,
% the bounds on measurement noise are minimized.
%
% The input data at time T (where T is the time horizon) is not
% used but the code requires the equal length for input and output, hence,
% if the input has length T-1, then it is allowed to add any value as the
% Tth element.
%
% Author: MI4Hybrid
% Date: Mar 29th, 2016

% Check if the system model is valid for this function.
if(strcmp(sys.mark,'poly')~=1)
    error('The system model must be a (non-switched) polynomial model.');
end

% Check if the uncertainty terms are specified and nonzeros.
if(isprop(sys,'d_coeffmat'))
    if(sum(sum(sys.d_coeffmat))==0)
        case1='cer';
    else
        case1='unc';
    end
else
    case1='cer';
end
% Obtain the system information.
n=size(sys.coeffmat,1); % state/output dimension
n_i=size(sys.degmat,2)-n; % input dimension
n_m=size(sys.coeffmat,2); % number of monomials

% Set up default values if some arguments are not specified.
num_arg=nargin;
if(num_arg==5)
    LS=zeros(n,1)-inf;
    US=zeros(n,1)+inf;
    param.scalingSW=0;
    param.relaxOrder=[];
    param.eqTolerance=[];
    verbose=1;
elseif(num_arg==7)
    param.scalingSW=0;
    param.relaxOrder=[];
    param.eqTolerance=[];
    verbose=1;
elseif(num_arg==8)
    verbose=1;
end

% Check if the argument "param" is valid.
if(ischar(param))
   param=[];
   param.scalingSW=0;
   param.relaxOrder=[];
   param.eqTolerance=[];
   warning(['The polynomial invalidation problem may be nontrivial'...
            ', the solver "SparsePOP" will be used, and the parameter '...
            'settings are changed.']);
end

if(isempty(verbose))
    verbose=1;
end
if(~isfield(param,{'scalingSW'}))
    param.scalingSW=0;
end
if(~isfield(param,{'relaxOrder'}))
    param.relaxOrder=[];
end
if(~isfield(param,{'eqTolerance'}))
    param.eqTolerance=10^-10;
end

% Set up default values for empty arguments.
if(isempty(BP))
    BP=zeros(n,1);
end
if(isempty(BM))
    BM=zeros(n,1);
end
if(isempty(LS))
    LS=zeros(n,1)-inf;
end
if(isempty(US))
    US=zeros(n,1)+inf;
end

% Convert scalars to vectors.
if(length(BP)==1&&n>1)
    BP=ones(n,1)*BP;
    warning(['Bound for process noise is a scalar, converted to a '...
        'vector with identical entries.']);
end
if(length(BM)==1&&n>1)
    BM=ones(n,1)*BM;
    warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
end
if(length(LS)==1&&n>1)
    LS=ones(n,1)*LS;
    warning(['Lower state bound is a scalar, converted to a vector '...
        'with identical entries.']);
end
if(length(US)==1&&n>1)
    US=ones(n,1)*US;
    warning(['Upper state bound is a scalar, converted to a vector '...
        'with identical entries.']);
end

% Check the bounds.
if(length(BP)~=n||~isvector(BP))
    error('The bound dimension for process noise is not correct.');
end
if(length(BM)~=n||~isvector(BM))
    error('The bound dimension for measurement noise is not correct.');
end
if(length(LS)~=n||~isvector(LS))
    error('The lower bound dimension for states is not correct.');
end
if(length(US)~=n||~isvector(US))
    error('The upper bound dimension for states is not correct.');
end
if(sum(BP<0)~=0)
    error('The process noise bound should not have negative value.');
end
if(sum(BM<0)~=0)
    error('The measurement noise bound should not have negative value.');
end

% Check if the input and output are consistent.
if(size(u,2)~=size(y,2))
    error('The lengths of input and output sequences are not consistent.');
end
% Check if the input is consistent with the model.
if(size(u,1)~=n_i)
    error('The input sequence is not consistent with the model.');
end
% Check if the output is consistent with the model.
if(size(y,1)~=n)
    error('The output sequence is not consistent with the model.');
end

% There are 2 different cases in total.
if(strcmp(case1,'cer'))
    mycase=1;
elseif(strcmp(case1,'unc'))
    mycase=2;
end

% Apply different functions for different cases.
switch mycase
    case 1
        [miflag,EPS]=...
            inv_poly_cer(sys,u,y,BP,BM,LS,US,param,verbose);
    case 2
        [miflag,EPS]=...
            inv_poly_unc(sys,u,y,BP,BM,LS,US,param,verbose);
end

end

%% Polynomial model invalidation without uncertainty.

function [miflag,EPS]=inv_poly_cer(sys,u,y,BP,BM,LS,US,param,verbose)

% Obtain the system information.
n=size(sys.coeffmat,1); % state/output dimension
n_i=size(sys.degmat,2)-n; % input dimension
n_m=size(sys.coeffmat,2); % number of monomials

% Obtain the time horizon.
T=size(y,2);

% Find normalization diagonal matrix for measurement noise bounds.
T_c = zeros(n,n);
for i = 1:n
    if BM(i)~=0
        T_c(i,i) = 1/BM(i);
    else
        T_c(i,i)=inf;
    end
end

% Define noise variables for the invalidation problem.
etap=sdpvar(n,T-1);
etam=sdpvar(n,T);
eps = sdpvar(1,1);
obj = eps;

% setting up constraints:
const=[];
for i=1:T-1
    current_var=y(:,i)-sys.Em*etam(:,i);
    for j=1:n_m % Calculate all monomials at time point i.
        temp=1;
        for l=1:n
            temp=temp*current_var(l)^sys.degmat(j,l);
        end
        for l=1:n_i
            temp=temp*u(l,i)^sys.degmat(j,l+n);
        end
        mono(j)=temp;
    end
    for j=1:n % Calculate polynomials using those monomials.
        poly(j,1)=sys.coeffmat(j,:)*mono';
    end
    % Set up constraints.
    const=[const,poly+sys.Ep*etap(:,i)-y(:,i+1)+sys.Em*etam(:,i+1)==0];
end

% measurement noise bound and state bound constraints.
for i=1:T
    const = [const, y(:,i)-etam(:,i)<=US];
    const = [const, y(:,i)-etam(:,i)>=LS];
end
for i = 1:T
    const=[const,norm(T_c*etam(:,i),inf)<=eps];
end

% process noise bound constraints.
for i=1:T-1
    for j=1:n
        const=[const,etap(j,i)<=BP(j),etap(j,i)>=-BP(j)];
    end
end

% Choose a solver and set up parameters.
ops=sdpsettings('solver','sparsepop','sparsepop.relaxOrder',param.relaxOrder,...
    'sparsepop.scalingSW',param.scalingSW,'sparsepop.eqTolerance', ... 
    param.eqTolerance,'verbose',verbose);

% Check the feasibility problem.
optimize(const,obj,ops);
EPS = double(eps);
% check if the model is invalidated or not
if EPS <= 1
    miflag = 0;  % if not invalidated
else
    miflag = 1;  % if invalidated
end
end

%% Polynomial model invalidation with uncertainty.

function [miflag,EPS]=inv_poly_unc(sys,u,y,BP,BM,LS,US,param,verbose)

% Obtain the system information.
n=size(sys.coeffmat,1); % state/output dimension
n_i=size(sys.degmat,2)-n; % input dimension
n_m=size(sys.coeffmat,2); % number of monomials

% no negative bounds on uncertainty
if(sum(sum(sys.d_coeffmat<0))~=0)
    error('The uncertainty bound should not have negative value.');
end

% Obtain the time horizon.
T=size(y,2);
% Find normalization diagonal matrix for measurement noise bounds.
T_c = zeros(n,n);
for i = 1:n
    if BM(i)~=0
        T_c(i,i) = 1/BM(i);
    else
        T_c(i,i)=inf;
    end
end
% Define noise variables for the invalidation problem.
etap=sdpvar(n,T-1);
etam=sdpvar(n,T);
eps = sdpvar(1,1);
uncer=sdpvar(n,n_m);
uncer(sys.d_coeffmat==0)=0;
obj = eps;

% setting up constraints:
const=[];
for i=1:T-1
    current_var=y(:,i)-sys.Em*etam(:,i);
    for j=1:n_m % Calculate all monomials at time point i.
        temp=1;
        for l=1:n
            temp=temp*current_var(l)^sys.degmat(j,l);
        end
        for l=1:n_i
            temp=temp*u(l,i)^sys.degmat(j,l+n);
        end
        mono(j)=temp;
    end
    for j=1:n % Calculate polynomials using those monomials.
        poly(j,1)=sys.coeffmat(j,:)*mono'+uncer(j,:)*mono';
    end
    % Set up constraints
    const=[const,poly+sys.Ep*etap(:,i)-y(:,i+1)+sys.Em*etam(:,i+1)==0];
end

% measurement noise bound and state bound constraints.
for i=1:T
        const = [const, y(:,i)-etam(:,i)<=US];
        const = [const, y(:,i)-etam(:,i)>=LS];
end
for i = 1:T
    const=[const,norm(T_c*etam(:,i),inf)<=eps];
end

% process noise bound constraints.
for i=1:T-1
    for j=1:n
        const=[const,etap(j,i)<=BP(j),etap(j,i)>=-BP(j)];
    end
end

% uncertainty bound constraints.
for i=1:n
    for j=1:n_m
        const=[const,uncer(i,j)<=sys.d_coeffmat(i,j),...
                uncer(i,j)>=-sys.d_coeffmat(i,j)];
    end
end

% Choose a solver and set up parameters.
ops=sdpsettings('solver','sparsepop','sparsepop.relaxOrder',param.relaxOrder,...
    'sparsepop.scalingSW',param.scalingSW,'sparsepop.eqTolerance', ... 
    param.eqTolerance,'verbose',verbose);

% Check the feasibility problem.
optimize(const,obj,ops);
EPS = double(eps);
% check if the model is invalidated or not
if EPS <= 1
    miflag = 0;  % if not invalidated
else
    miflag = 1;  % if invalidated
end
end


