function [miflag,primal,dual]=invalidation_poly(sys,input,output,...
    pn_bound,mn_bound,state_bound,param)

% This function will test if a non-switched polynomial model is invalidated
% or not based on the input, output and noise bounds. The function requires
% the SparsePOP package and SDP solvers such as SEDUMI.
%
% Arguments:
%   sys -- a user defined polynomial model (see PolyModel.m)
%   input -- the input to the system, an n_i-by-T matrix where n_i is the
%            input dimension and T is the time horizon
%   output -- the output from the system, an n_y-by-T matrix where n_y is
%             the output dimension and T is the time horizon
%   pn_bound -- n-D vector specifying the infinity norm bounds for process
%               noise where n is the state dimension
%   mn_bound -- n_y-D vector specifying the infinity norm bounds for
%               measurement noise where n_y is the output dimension
%   state_bound -- n-D vector specifying the infinity norm bounds for
%                  states where n is the state dimension
%   param -- parameter structure for SparsePOP
%            User can refine the solution obtained from the SDP relaxation
%            by setting the parameter:
%               param.POPsolver = 'active-set';
%               param.POPsolver = 'interior-point';
%               param.POPsolver = 'trust-region-reflective';
%               param.POPsolver = 'lsqnonlin';
%            See SparsePOP userGuide.pdf for more information about
%               param.relaxOrder
%
% Outputs:
%   miflag -- a flag indicating if the model is invalidated (0 for not
%             invalidated and 1 for invalidated)
%   primal -- solution to the optimization
%   dual -- the dual corresponding to the primal: model is not invalidated
%           if nagative and model is invalidated if positive
%
% Note: n_y is equivalent to n because the polynomial model assumes that
% the output is the measurement of states.
%
% Syntax:
%   [miflag,primal,dual]=invalidation_poly(sys,input,output,pn_bound,...
%                           mn_bound);
%   [miflag,primal,dual]=invalidation_poly(sys,input,output,pn_bound,...
%                           mn_bound,state_bound);
%   [miflag,primal,dual]=invalidation_poly(sys,input,output,pn_bound,...
%                           mn_bound,state_bound,param);
%
% Author: MI4Hybrid
% Date: July 24th, 2015

% Check if the system model is valid for this function.
if(strcmp(sys.mark,'poly')~=1)
    error('The system model must be a (non-switched) polynomial model.');
end

% Obtain system parameter.
n=size(sys.coeffmat,1); % state/output dimension
n_i=size(sys.degmat,2)-n; % input dimension

% Set up default values if some arguments are not specified.
num_arg=nargin;
if(num_arg==5)
    state_bound=zeros(n,1)+inf;
    param.eqTolerance=10^-10;
    param.scalingSW=1;
elseif(num_arg==6)
    param.eqTolerance=10^-10;
    param.scalingSW=1;
end

% Set up default values for empty arguments.
if(isempty(pn_bound))
    pn_bound=zeros(n,1);
end
if(isempty(mn_bound))
    mn_bound=zeros(n,1);
end
if(isempty(state_bound))
    state_bound=zeros(n,1)+inf;
end
if(isempty(param))
    param.eqTolerance=10^-10;
    param.scalingSW=1;
end

% Convert scalars to vectors.
if(length(pn_bound)==1&&n>1)
    pn_bound=ones(n,1)*pn_bound;
    warning(['Bound for process noise is a scalar, converted to a '...
        'vector with identical entries.']);
end
if(length(mn_bound)==1&&n>1)
    mn_bound=ones(n,1)*mn_bound;
    warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
end
if(length(state_bound)==1&&n>1)
    state_bound=ones(n,1)*state_bound;
    warning(['State bound is a scalar, converted to a vector with '...
        'identical entries.']);
end

% Check the bounds.
if(length(pn_bound)~=n||~isvector(pn_bound))
    error('The number of bounds for process noise is not correct.');
end
if(length(mn_bound)~=n||~isvector(mn_bound))
    error('The number of bounds for measurement noise is not correct.');
end
if(length(state_bound)~=n||~isvector(state_bound))
    error('The number of bounds for states is not correct.');
end

% Check if the input and output are consistent.
if(size(input,2)~=size(output,2))
    error('The input length and output length are not consistent.');
end
% Check if the input is consistent with the model.
if(size(input,1)~=n_i)
    error('The input is not consistent with the model.');
end
% Check if the output is consistent with the model.
if(size(output,1)~=n)
    error('The output is not consistent with the model.');
end

% Obtain the time horizon after checking.
T=size(input,2); % time horizon

% Construct the symbolic polynomial functions for the system model.
PolyVar=mpvar('v',[n_i+n 1]);
P=polynomial(sys.coeffmat',sys.degmat,PolyVar.varname,[n,1]);

% Construct polynomials with I/O data and noise variables.
noise_var1=mpvar('m',[T*n 1]);
noise_var2=mpvar('p',[(T-1)*n 1]);
new_P=mpvar('u',[n T-1]);
for i=1:T-1
    var_buffer=output(:,i)-sys.Em*noise_var1((i-1)*n+1:i*n);
    new_P(:,i)=subs(P,PolyVar,[var_buffer;input(:,i)]);
    new_P(:,i)=new_P(:,i)-output(:,i+1)+sys.Ep*noise_var2((i-1)*n+1:i*n)...
        +sys.Em*noise_var1(i*n+1:(i+1)*n);
end

% Construct the objective function (sum of squares polynomial).
obj=sum(sum(new_P.^2));

% The degree matrix and coefficient matirx of the objective function.
coeffmat=full(obj.coefficient);
coeffmat=coeffmat/max(abs(coeffmat));
degmat=full(obj.degmat);

% Construct the SparsePOP data format for the objective function.
objPoly.typeCone=1;
objPoly.dimVar=length(obj.varname);
objPoly.degree=max(sum(degmat,2));
objPoly.noTerms=length(coeffmat);
objPoly.supports=degmat;
objPoly.coef=coeffmat;

% Obtain variable information of the objective function.
variables=obj.varname;
len=length(variables);

% Check the norm types for noise.
for i=1:n
    if(sys.mn_norm(i)~=inf)
        warning(['At least one dimension of measurement noise is '...
            'specified with a norm type which is not the infinity norm.'...
            ' The function will automatically apply the infinity norm.']);
        break
    elseif(sys.pn_norm(i)~=inf)
        warning(['At least one dimension of process noise is specified'...
            ' with a norm type which is not the infinity norm. The '...
            'function will automatically apply the infinity norm.']);
        break
    end
end

% Construct the SparsePOP data format for noise constraints.
ineqPolySys=[];
lbd=zeros(1,len); % pre-allocate memory
ubd=zeros(1,len); % pre-allocate memory
for i=1:len
    if(strcmp(variables{i}(1),'m'))
        suffix=str2double(variables{i}(3:end));
        idx=rem(suffix,n);
        time=floor(suffix/n);
        if(idx==0)
            idx=n;
        end
        lbd_candidate=-state_bound(idx)+output(idx,time);
        if(lbd_candidate<=-mn_bound(idx))
            lbd(i)=-mn_bound(idx);
        else
            lbd(i)=lbd_candidate;
        end
        ubd_candidate=state_bound(idx)+output(idx,time);
        if(ubd_candidate>=mn_bound(idx))
            ubd(i)=mn_bound(idx);
        else
            ubd(i)=ubd_candidate;
        end
    elseif(strcmp(variables{i}(1),'p'))
        idx=rem(str2double(variables{i}(3:end)),n);
        if(idx==0)
            idx=n;
        end
        lbd(i)=-pn_bound(idx);
        ubd(i)=pn_bound(idx);
    end
end

% Apply the SparsePOP.
[~,SDPobjValue,POP,~,~,~]=sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
primal=SDPobjValue;
dual=POP.objValue;
if(primal>0)
    miflag=1;
else
    miflag=0;
end

end