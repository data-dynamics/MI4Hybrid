%% The main function for T-detectability.

function [tdflag,epsilon]=tdetect_poly(sys,sysf,T,BM,BMf,LS,LSf,US,USf,...
    LU,LUf,UU,UUf,BP,BPf,param,verbose)

% This function uses noise, state, input bounds (and uncertainty bounds, if
% any) to determine if a fault of a (non-switched) polynomial model is
% detectable in T steps. The returned value can also indicate how "far" the
% fault is from being not detectable.
%
% Input:
%   sys -- a user defined polynomial model (a priori)
%   sysf -- a user defined fault polynomial model
%             (see PolyModel.m for details about polynomial models)
%   T -- time horizon (for T-detectability)
%   BP -- n-D vector specifying the process noise bound of "sys"
%   BM -- n-D vector specifying the measurement noise bound of "sys"
%   BPf -- n-D vector specifying the process noise bound of "sysf"
%   BMf -- n-D vector specifying the measurement noise bound of "sysf"
%   LS -- n-D vector specifying the lower bound of states of "sys"
%   US -- n-D vector specifying the upper bound of states of "sys"
%   LSf -- n-D vector specifying the lower bound of states of "sysf"
%   USf -- n-D vector specifying the upper bound of states of "sysf"
%             (n is the state/output dimension)
%   LU -- ni-D vector specifying the lower bound of inputs of "sys"
%   UU -- ni-D vector specifying the upper bound of inputs of "sys"
%   LUf -- ni-D vector specifying the lower bound of inputs of "sysf"
%   UUf -- ni-D vector specifying the upper bound of inputs of "sysf"
%             (ni is the input dimension)
%   param -- parameter structure for the solver "SparsePOP"
%            User can refine the solution obtained from the SDP relaxation
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
%   tdflag -- a flag indicating if the fault is T-detectable (0 for not
%             T-detectable and 1 for T-detectable)
%   epsilon -- the optimal value for the objective function. If this value
%              is NaN, it means that the optimization problem is infeasible
%              (and the model is invalidated).
%
% Syntax:
%  [tdflag,epsilon]=tdet_poly(sys,sysf,T,BM,BMf);
%  [tdflag,epsilon]=tdet_poly(sys,sysf,T,BM,BMf,LS,LSf,US,USf);
%  [tdflag,epsilon]=tdet_poly(sys,sysf,T,BM,BMf,LS,LSf,US,USf,LU,LUf,UU,...
%                              UUf);
%  [tdflag,epsilon]=tdet_poly(sys,sysf,T,BM,BMf,LS,LSf,US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf);
%  [tdflag,epsilon]=tdet_poly(sys,sysf,T,BM,BMf,LS,LSf,US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf,param);
%  [tdflag,epsilon]=tdet_poly(sys,sysf,T,BM,BMf,LS,LSf,US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf,param,verbose);
%
% Note: the objective function will be different for different cases (e.g.
% with or without uncertainty). Details are discussed in the paper "Model
% (In)validation and Fault Detection for Systems with Polynomial
% State-Space Models" by F. Harirchi, Z. Luo, and N. Ozay.
%
% Author: MI4Hybrid
% Date: Oct 27th, 2015

% Check if the system model is valid for this function.
if(strcmp(sys.mark,'poly')~=1)
    error('The a priori model must be a (non-switched) polynomial model.');
end
if(strcmp(sysf.mark,'poly')~=1)
    error('The fault model must be a (non-switched) polynomial model.');
end

% what about the uncertainty for the fault model?
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

m_degree = max(sum(sys.degmat,2));
% Set up default values if some arguments are not specified.
num_arg=nargin;
if(num_arg==5)
    LS=zeros(n,1)-inf;
    LSf=zeros(n,1)-inf;
    US=zeros(n,1)+inf;
    USf=zeros(n,1)+inf;
    LU=zeros(n_i,1)-inf;
    LUf=zeros(n_i,1)-inf;
    UU=zeros(n_i,1)+inf;
    UUf=zeros(n_i,1)+inf;
    BP=zeros(n,1);
    BPf=zeros(n,1);
    param.scalingSW=0;
    param.relaxOrder=[];
    param.eqTolerance=[];
    param.POPsolver=[];
    verbose=1;
elseif(num_arg==9)
    LU=zeros(n_i,1)-inf;
    LUf=zeros(n_i,1)-inf;
    UU=zeros(n_i,1)+inf;
    UUf=zeros(n_i,1)+inf;
    BP=zeros(n,1);
    BPf=zeros(n,1);
    param.scalingSW=0;
    param.relaxOrder=[];
    param.eqTolerance=1e-10;
    param.POPsolver=[];
    verbose=1;
elseif(num_arg==13)
    BP=zeros(n,1);
    BPf=zeros(n,1);
    param.scalingSW=0;
    param.relaxOrder=floor(m_degree/2)+1;
    param.eqTolerance=1e-10;
    param.POPsolver=[];
    verbose=1;
elseif(num_arg==15)
    param.scalingSW=0;
    param.relaxOrder=floor(m_degree/2)+1;
    param.eqTolerance=1e-10;
    param.POPsolver=[];
    verbose=1;
elseif(num_arg==16)
    verbose=1;
end

% Check if the system model is valid for this function.
if(strcmp(sys.mark,'poly')~=1)
    error('The a priori model must be a (non-switched) polynomial model.');
end
if(strcmp(sysf.mark,'poly')~=1)
    error('The fault model must be a (non-switched) polynomial model.');
end

% Obtain system parameters.
n=size(sys.coeffmat,1); % state/output dimension
if(n~=size(sysf.coeffmat,1))
    error('The two models must have the same state/output dimension.');
end
n_i=size(sys.degmat,2)-n; % input dimension
if(n_i~=size(sysf.degmat,2)-n)
    error('The two models must have the same input dimension.');
end
n_m=size(sys.coeffmat,2); % number of monomials of "sys"
n_mf=size(sysf.coeffmat,2); % number of monomials of "sysf"

% Set up default values if some arguments are not specified.
num_arg=nargin;
if(num_arg==5)
    LS=zeros(n,1)-inf;
    LSf=zeros(n,1)-inf;
    US=zeros(n,1)+inf;
    USf=zeros(n,1)+inf;
    LU=zeros(n_i,1)-inf;
    LUf=zeros(n_i,1)-inf;
    UU=zeros(n_i,1)+inf;
    UUf=zeros(n_i,1)+inf;
    BP=zeros(n,1);
    BPf=zeros(n,1);
    param.scalingSW=0;
    param.relaxOrder=floor(m_degree/2)+1;
    param.eqTolerance=1e-10;
    verbose=1;
elseif(num_arg==9)
    LU=zeros(n_i,1)-inf;
    LUf=zeros(n_i,1)-inf;
    UU=zeros(n_i,1)+inf;
    UUf=zeros(n_i,1)+inf;
    BP=zeros(n,1);
    BPf=zeros(n,1);
    param.scalingSW=0;
    param.relaxOrder=floor(m_degree/2)+1;
    param.eqTolerance=1e-10;
    verbose=1;
elseif(num_arg==13)
    BP=zeros(n,1);
    BPf=zeros(n,1);
    param.scalingSW=0;
    param.relaxOrder=floor(m_degree/2)+1;
    param.eqTolerance=1e-10;
    verbose=1;
elseif(num_arg==15)
    param.scalingSW=0;
    param.relaxOrder=floor(m_degree/2)+1;
    param.eqTolerance=1e-10;
    verbose=1;
elseif(num_arg==16)
    verbose=1;
end

% Set up default values for empty arguments.
if(isempty(LS))
    LS=zeros(n,1)-inf;
end
if(isempty(LSf))
    LSf=zeros(n,1)-inf;
end
if(isempty(US))
    US=zeros(n,1)+inf;
end
if(isempty(USf))
    USf=zeros(n,1)+inf;
end
if(isempty(LU))
    LU=zeros(n,1)-inf;
end
if(isempty(LUf))
    LUf=zeros(n,1)-inf;
end
if(isempty(UU))
    UU=zeros(n,1)+inf;
end
if(isempty(UUf))
    UUf=zeros(n,1)+inf;
end
if(isempty(BP))
    BP=zeros(n,1);
end
if(isempty(BPf))
    BPf=zeros(n,1);
end
if(isempty(param))
    param.scalingSW=0;
    param.relaxOrder=floor(m_degree/2)+1;
    param.eqTolerance=1e-10;
    param.POPsolver=[];
end
if(isempty(verbose))
    verbose=1;
end
if(~isfield(param,{'scalingSW'}))
    param.scalingSW=0;
end
if(~isfield(param,{'relaxOrder'}))
    param.relaxOrder=floor(m_degree/2)+1;
end
if(~isfield(param,{'eqTolerance'}))
    param.eqTolerance=1e-10;
end

% Convert scalars to vectors.
if(length(BP)==1&&n>1)
    BP=ones(n,1)*BP;
    warning(['Bound for process noise of the a priori model is a scalar'...
        ', converted to a vector with identical entries.']);
end
if(length(BPf)==1&&n>1)
    BPf=ones(n,1)*BPf;
    warning(['Bound for process noise of the fault model is a scalar, '...
        'converted to a vector with identical entries.']);
end
if(length(BM)==1&&n>1)
    BM=ones(n,1)*BM;
    warning(['Bound for measurement noise of the a priori model is a '...
        'scalar, converted to a vector with identical entries.']);
end
if(length(BMf)==1&&n>1)
    BMf=ones(n,1)*BMf;
    warning(['Bound for measurement noise of the fault model is a '...
        'scalar, converted to a vector with identical entries.']);
end
if(length(LS)==1&&n>1)
    LS=ones(n,1)*LS;
    warning(['Lower state bound of the a priori model is a scalar, '...
        'converted to a vector with identical entries.']);
end
if(length(LSf)==1&&n>1)
    LSf=ones(n,1)*LSf;
    warning(['Lower state bound of the fault model is a scalar, '...
        'converted to a vector with identical entries.']);
end
if(length(US)==1&&n>1)
    US=ones(n,1)*US;
    warning(['Upper state bound of the a priori model is a scalar, '...
        'converted to a vector with identical entries.']);
end
if(length(USf)==1&&n>1)
    USf=ones(n,1)*USf;
    warning(['Upper state bound of the fault model is a scalar, '...
        'converted to a vector with identical entries.']);
end
if(length(LU)==1&&n_i>1)
    LU=ones(n_i,1)*LU;
    warning(['Lower input bound of the a priori model is a scalar, '...
        'converted to a vector with identical entries.']);
end
if(length(LUf)==1&&n_i>1)
    LUf=ones(n_i,1)*LUf;
    warning(['Lower input bound of the fault model is a scalar, '...
        'converted to a vector with identical entries.']);
end
if(length(UU)==1&&n_i>1)
    UU=ones(n_i,1)*UU;
    warning(['Upper input bound of the a priori model is a scalar, '...
        'converted to a vector with identical entries.']);
end
if(length(UUf)==1&&n_i>1)
    UUf=ones(n_i,1)*UUf;
    warning(['Upper input bound of the fault model is a scalar, '...
        'converted to a vector with identical entries.']);
end

% Check the bounds.
if(length(BP)~=n||~isvector(BP))
    error(['The bound dimension for process noise of the a priori model'...
        ' is not correct.']);
end
if(length(BPf)~=n||~isvector(BPf))
    error(['The bound dimension for process noise of the fault model '...
        'is not correct.']);
end
if(length(BM)~=n||~isvector(BM))
    error(['The bound dimension for measurement noise of the a priori '...
        'model is not correct.']);
end
if(length(BMf)~=n||~isvector(BMf))
    error(['The bound dimension for measurement noise of the fault '...
        'model is not correct.']);
end
if(length(LS)~=n||~isvector(LS))
    error(['The lower bound dimension for states of the a priori model'...
        ' is not correct.']);
end
if(length(LSf)~=n||~isvector(LSf))
    error(['The lower bound dimension for states of the fault model is'...
        ' not correct.']);
end
if(length(US)~=n||~isvector(US))
    error(['The upper bound dimension for states of the a priori model'...
        ' is not correct.']);
end
if(length(USf)~=n||~isvector(USf))
    error(['The upper bound dimension for states of the fault model is'...
        ' not correct.']);
end
if(length(LU)~=n||~isvector(LU))
    error(['The lower bound dimension for inputs of the a priori model'...
        ' is not correct.']);
end
if(length(LUf)~=n||~isvector(LUf))
    error(['The lower bound dimension for inputs of the fault model is'...
        ' not correct.']);
end
if(length(UU)~=n||~isvector(UU))
    error(['The upper bound dimension for inputs of the a priori model'...
        ' is not correct.']);
end
if(length(UUf)~=n||~isvector(UUf))
    error(['The upper bound dimension for inputs of the fault model is'...
        ' not correct.']);
end
if(sum(BP<0)~=0)
    error(['The process noise bound of the a priori model should not '...
        'have negative value.']);
end
if(sum(BPf<0)~=0)
    error(['The process noise bound of the fault model should not have '...
        'negative value.']);
end
if(sum(BM<0)~=0)
    error(['The measurement noise bound of the a priori model should '...
        'not have negative value.']);
end
if(sum(BMf<0)~=0)
    error(['The measurement noise bound of the fault model should not '...
        'have negative value.']);
end


% Apply different functions for different cases.
switch case1
    case 'cer'
        [tdflag,epsilon]=tdet_poly_cer(sys,sysf,T,BM,BMf,LS,...
            LSf,US,USf,LU,LUf,UU,UUf,BP,BPf,param,verbose);
    case 'unc'
        [tdflag,epsilon]=tdet_poly_unc(sys,sysf,T,BM,BMf,LS,...
            LSf,US,USf,LU,LUf,UU,UUf,BP,BPf,param,verbose);
end

end

%% T-detectability without uncertainty.

function [tdflag,epsilon]=tdet_poly_cer(sys,sysf,T,BM,BMf,LS,...
    LSf,US,USf,LU,LUf,UU,UUf,BP,BPf,param,verbose)

% This function uses noise, state, and input bounds to determine if a fault
% of a non-switched polynomial model is detectable in T steps. The returned
% value can also indicate how "far" the fault is from being not detectable.
%
% Input:
%   sys -- a user defined polynomial model (a priori)
%   sysf -- a user defined fault polynomial model
%             (see PolyModel.m for details about polynomial models)
%   T -- time horizon (for T-detectability)
%   BP -- n-D vector specifying the process noise bound of "sys"
%   BM -- n-D vector specifying the measurement noise bound of "sys"
%   BPf -- n-D vector specifying the process noise bound of "sysf"
%   BMf -- n-D vector specifying the measurement noise bound of "sysf"
%   LS -- n-D vector specifying the lower bound of states of "sys"
%   US -- n-D vector specifying the upper bound of states of "sys"
%   LSf -- n-D vector specifying the lower bound of states of "sysf"
%   USf -- n-D vector specifying the upper bound of states of "sysf"
%             (n is the state/output dimension)
%   LU -- ni-D vector specifying the lower bound of inputs of "sys"
%   UU -- ni-D vector specifying the upper bound of inputs of "sys"
%   LUf -- ni-D vector specifying the lower bound of inputs of "sysf"
%   UUf -- ni-D vector specifying the upper bound of inputs of "sysf"
%             (ni is the input dimension)
%   param -- parameter structure for the solver "SparsePOP"
%            User can refine the solution obtained from the SDP relaxation
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
%   tdflag -- a flag indicating if the fault is T-detectable (0 for not
%             T-detectable and 1 for T-detectable)
%   epsilon -- optimal solution for a dimension of the measurement noise
%            bound of "sys" (can be multiple dimensions if those dimensions
%            of noise are bounded by the same value). If this value is NaN,
%            it means that the optimization problem is not fesible (and of
%            course, the model is invalidated).
%
% Syntax:
%   [tdflag,epsilon]=tdet_poly_cer(sys,sysf,T,BM,BMf);
%   [tdflag,epsilon]=tdet_poly_cer(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf);
%   [tdflag,epsilon]=tdet_poly_cer(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf,LU,LUf,UU,...
%                              UUf);
%   [tdflag,epsilon]=tdet_poly_cer(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf);
%   [tdflag,epsilon]=tdet_poly_cer(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf,param);
%   [tdflag,epsilon]=tdet_poly_cer(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf,param,verbose);
%
%
% Author: MI4Hybrid
% Date: Oct 26th, 2015



% Find repeated values in the vector of measurement noise bound.
% Find normalization diagonal matrix for measurement noise bounds.
n=size(sys.coeffmat,1); % state/output dimension
n_i=size(sys.degmat,2)-n; % input dimension
n_m=size(sys.coeffmat,2); % number of monomials of "sys"
n_mf=size(sysf.coeffmat,2); % number of monomials of "sysf"
T_c = zeros(n,n);
for i = 1:n
    if BM(i)~=0
        T_c(i,i) = 1/BM(i);
    else
        T_c(i,i)=inf;
    end
end
% Define the objective function.
eps=sdpvar(1,1);
obj=eps;

% Define measurement noise variables.
etam=sdpvar(n,T);
etam(BM==0,:)=0;
etamb=sdpvar(n,T);
etamb(BMf==0,:)=0;

% Define process noise variables.
etap=sdpvar(n,T-1);
etap(BP==0,:)=0;
etapb=sdpvar(n,T-1);
etapb(BPf==0,:)=0;

% Define input/output variables.
ubH=min([UU UUf],[],2); % upper input bound
ubL=max([LU LUf],[],2); % lower input bound
temp=abs(ubH)+abs(ubL);
u=sdpvar(n_i,T-1);
u(temp==0,:)=0;
y=sdpvar(n,T);


% Define the constraints using model equations of "sys".
const=[];
sdpvar mono(1,2)
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


% Define the constraints using model equations of "sysf".
sdpvar monof(1,2)
for i=1:T-1
    current_var=y(:,i)-sysf.Em*etamb(:,i);
    for j=1:n_mf % Calculate all monomials at time point i.
        temp=1;
        for l=1:n
            temp=temp*current_var(l)^sysf.degmat(j,l);
        end
        for l=1:n_i
            temp=temp*u(l,i)^sysf.degmat(j,l+n);
        end
        monof(j)=temp;
    end
    for j=1:n % Calculate polynomials using those monomials.
        polyf(j,1)=sysf.coeffmat(j,:)*monof';
    end
    % Set up constraints.
    const=[const,polyf+sysf.Ep*etapb(:,i)-y(:,i+1)+sysf.Em*etamb(:,i+1)==0];
end

% Define the constraints using noise and state bounds.
for i=1:T
    const=[const,norm(T_c*etam(:,i),inf)<=eps];
    for j=1:n
        const=[const,y(j,i)-etam(j,i)<=US(j),...
            y(j,i)-etam(j,i)>=LS(j),y(j,i)-etamb(j,i)<=USf(j),...
            y(j,i)-etamb(j,i)>=LSf(j),etamb(j,i)<=BMf(j),...
            etamb(j,i)>=-BMf(j)];
    end
end

for i=1:T-1
    for j=1:n
        const=[const,etap(j,i)<=BP(j),etap(j,i)>=-BP(j),...
            etapb(j,i)<=BPf(j),etapb(j,i)>=-BPf(j)];
    end
end

% Define the constraints using input bound.
for i=1:T-1
    for j=1:n_i
        const=[const,u(j,i)<=ubH(j),u(j,i)>=ubL(j,i)];
    end
end



% Set up the solver settings, and solve the optimization problem.
ops=sdpsettings('verbose',verbose,'solver','sparsepop',...
    'sparsepop.relaxOrder',3,'sparsepop.eqTolerance',param.eqTolerance);

optimize(const,obj,ops);
epsilon=double(eps);
if(epsilon>1||isnan(epsilon))
    tdflag=1;
else
    tdflag=0;
end

end

%% T-detectability with uncertainty.
function [tdflag,epsilon]=tdet_poly_unc(sys,sysf,T,BM,BMf,LS,...
    LSf,US,USf,LU,LUf,UU,UUf,BP,BPf,param,verbose)

% This function uses noise, state, and input bounds to determine if a fault
% of a non-switched polynomial model (with uncertainty) is detectable in T
% steps. The returned value can also indicate how "far" the fault is from
% being not detectable.
%
% Input:
%   sys -- a user defined polynomial model (a priori)
%   sysf -- a user defined fault polynomial model
%             (see PolyModel.m for details about polynomial models)
%   T -- time horizon (for T-detectability)
%   BP -- n-D vector specifying the process noise bound of "sys"
%   BM -- n-D vector specifying the measurement noise bound of "sys"
%   BPf -- n-D vector specifying the process noise bound of "sysf"
%   BMf -- n-D vector specifying the measurement noise bound of "sysf"
%   LS -- n-D vector specifying the lower bound of states of "sys"
%   US -- n-D vector specifying the upper bound of states of "sys"
%   LSf -- n-D vector specifying the lower bound of states of "sysf"
%   USf -- n-D vector specifying the upper bound of states of "sysf"
%             (n is the state/output dimension)
%   LU -- ni-D vector specifying the lower bound of inputs of "sys"
%   UU -- ni-D vector specifying the upper bound of inputs of "sys"
%   LUf -- ni-D vector specifying the lower bound of inputs of "sysf"
%   UUf -- ni-D vector specifying the upper bound of inputs of "sysf"
%             (ni is the input dimension)
%   param -- parameter structure for the solver "SparsePOP"
%            User can refine the solution obtained from the SDP relaxation
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
%   tdflag -- a flag indicating if the fault is T-detectable (0 for not
%             T-detectable and 1 for T-detectable)
%   epsilon -- the optimal solution for an uncertainty term (for multiple
%            uncertainty terms if those terms are assumed to have the same
%            bound) of "sys". If this value is NaN, then the optimization
%            problem is infesible (and thus, the model is invalidated).
%
% Syntax:
%   [tdflag,epsilon]=tdet_poly_unc(sys,sysf,T,BM,BMf);
%   [tdflag,epsilon]=tdet_poly_unc(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf);
%   [tdflag,epsilon]=tdet_poly_unc(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf,LU,LUf,UU,...
%                              UUf);
%   [tdflag,epsilon]=tdet_poly_unc(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf);
%   [tdflag,epsilon]=tdet_poly_unc(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf,param);
%   [tdflag,epsilon]=tdet_poly_unc(sys,sysf,T,BM,BMf,LS,LSf,...
%                              US,USf,LU,LUf,UU,...
%                              UUf,BP,BPf,param,verbose);
%
%
% Author: MI4Hybrid
% Date: Oct 27th, 2015



% Obtain system parameters.
n=size(sys.coeffmat,1); % state/output dimension
n_i=size(sys.degmat,2)-n; % input dimension
n_m=size(sys.coeffmat,2); % number of monomials of "sys"
n_mf=size(sysf.coeffmat,2); % number of monomials of "sysf"

T_c = zeros(n,n);
for i = 1:n
    if BM(i)~=0
        T_c(i,i) = 1/BM(i);
    else
        T_c(i,i)=inf;
    end
end
% % Find repeated values in the uncertainty matrix of "sys".
% unc=unique(sys.d_coeffmat,'stable');
% unc_length=length(unc);
% count=0;
% for i=1:unc_length
%     rep=length(find(sys.d_coeffmat==unc(i)));
%     if(rep>count&&unc(i)~=0)
%         count=rep; % the number of entries for the critical value
%         value=unc(i); % the critical value
%     end
% end
% eps_idx=find(sys.d_coeffmat==value); % the corresponding indices

% Define the objective function.
eps=sdpvar(1,1);
obj=eps;

% Define measurement noise variables.
etam=sdpvar(n,T);
etam(BM==0,:)=0;
etamb=sdpvar(n,T);
etamb(BMf==0,:)=0;

% Define process noise variables.
etap=sdpvar(n,T-1);
etap(BP==0,:)=0;
etapb=sdpvar(n,T-1);
etapb(BPf==0,:)=0;

% Define input/output variables.
ubH=min([UU UUf],[],2); % upper input bound
ubL=max([LU LUf],[],2); % lower input bound
temp=abs(ubH)+abs(ubL);
u=sdpvar(n_i,T-1);
u(temp==0,:)=0;
y=sdpvar(n,T);

% Define the uncertainty variables.
uncer=sdpvar(n,n_m);
uncer(sys.d_coeffmat==0)=0;
uncerb=sdpvar(n,n_mf);
uncerb(sysf.d_coeffmat==0)=0;
sdpvar mono(1,n_m)
% Define the constraints using model equations of "sys".
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
    % Set up constraints.
    const=[const,poly+sys.Ep*etap(:,i)-y(:,i+1)+sys.Em*etam(:,i+1)==0];
end

% Define the constraints using model equations of "sysf".
sdpvar monof(1,n_m)
for i=1:T-1
    current_var=y(:,i)-sysf.Em*etamb(:,i);
    for j=1:n_mf % Calculate all monomials at time point i.
        temp=1;
        for l=1:n
            temp=temp*current_var(l)^sysf.degmat(j,l);
        end
        for l=1:n_i
            temp=temp*u(l,i)^sysf.degmat(j,l+n);
        end
        monof(j)=temp;
    end
    for j=1:n % Calculate polynomials using those monomials.
        polyf(j,1)=sysf.coeffmat(j,:)*monof'+uncerb(j,:)*monof';
    end
    % Set up constraints.
    const=[const,polyf+sysf.Ep*etapb(:,i)-y(:,i+1)+sysf.Em*etamb(:,i+1)==0];
end

% Define the constraints using noise and state bounds.
for i=1:T
    for j=1:n
        const=[const,y(j,i)-etam(j,i)<=US(j),...
            y(j,i)-etam(j,i)>=LS(j),y(j,i)-etamb(j,i)<=USf(j),...
            y(j,i)-etamb(j,i)>=LSf(j),etamb(j,i)<=BMf(j),...
            etamb(j,i)>=-BMf(j),etam(j,i)<=BM(j),etam(j,i)>=-BM(j)];
    end
end

for i=1:T-1
    for j=1:n
        const=[const,etap(j,i)<=BP(j),etap(j,i)>=-BP(j),...
            etapb(j,i)<=BPf(j),etapb(j,i)>=-BPf(j)];
    end
end

% Define the constraints using input bound.
for i=1:T-1
    for j=1:n_i
        const=[const,u(j,i)<=ubH(j),u(j,i)>=ubL(j,i)];
    end
end

% Define the constraints using uncertainty bounds.
for i=1:n
    for j=1:n_m
        const=[const,uncerb(i,j)<=sysf.d_coeffmat(i,j),...
            uncerb(i,j)>=-sysf.d_coeffmat(i,j)];
        const=[const,uncer(i,j)<=sys.d_coeffmat(i,j),...
                uncer(i,j)>=-sys.d_coeffmat(i,j)];
    end
end

% Define the constraints associated with the objective function.
for i = 1:T
    const=[const,norm(T_c*etam(:,i),inf)<=eps];
end

% Set up the solver settings, and solve the optimization problem.
ops=sdpsettings('verbose',verbose,'solver','sparsepop',...
    'sparsepop.relaxOrder',param.relaxOrder,'sparsepop.eqTolerance',param.eqTolerance);

optimize(const,obj,ops);
epsilon=double(eps);
if(epsilon>1||isnan(epsilon))
    tdflag=1;
else
    tdflag=0;
end

end