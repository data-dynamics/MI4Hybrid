function result=invalidation_ss(sys,input,output,pn_bound,mn_bound,...
    state_bound,solver)

% This function applies the invalidation algorithm for state-space models
% based on the system input and output. This invlidation algorithm is made
% to check the feasibility of a linear equation in the form of M*x=b where
% x is the variable vector constructed by noise and state variables.
%
% Input/parameters:
%   sys -- the system model
%   input -- the input sequence
%   output -- the output sequence
%   pn_bound -- bound for process noise
%   mn_bound -- bound for measurement noise
%   state_bound -- bound for states
%   solver -- solver to be used to solve the optimization problem
%           -- 'cplex': uses CPLEX solver from IBM (needs to be installed)
%               for more information on cplex see: "http://www-03.ibm.com
%               /software/products/en/ibmilogcpleoptistud"
%           -- 'mosek': use Mosek as the solver (needs to be installed)
%               for more information on Mosek see: "https://www.mosek.com/"
% Output:
%   result -- return "true" if the system is validated, otherwise "false"
%
% Syntax:
%   result=invalidation_ss(sys,input,output,pn_bound,mn_bound,state_bound);
%   result=Invalidation_ss(sys,input,output,pn_bound,mn_bound,...
%                           state_bound,solver);
%
% Author: MI4Hybrid
% Date: July 8th, 2015

% Check if the system model is valid for this function.
if(strcmp(sys.mark,'ss')~=1)
    error('The system model must be a non-switched state-space model.');
end

% Obtain system model information.
T=size(input,2); % time horizon
n=size(sys.mode.C,2); % state dimension
n_y=size(sys.mode.C,1); % output dimension

% Use the default solver if it is not specified.
if(nargin==6)
    solver='cplex';
end

% Set up default values for empty paramters.
if(isempty(pn_bound))
    pn_bound=zeros(n,1);
end
if(isempty(mn_bound))
    mn_bound=zeros(n_y,1);
end
if(isempty(state_bound))
    state_bound=zeros(n,1)+inf;
end
if(isempty(solver))
    solver='cplex';
end

% Convert scalars to vectors.
if(length(pn_bound)==1&&n>1)
    pn_bound=ones(n,1)*pn_bound;
    warning(['Bound for process noise is a scalar, converted to a '...
        'vector with identical entries.']);
end
if(length(mn_bound)==1&&n_y>1)
    mn_bound=ones(n_y,1)*mn_bound;
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
if(length(mn_bound)~=n_y||~isvector(mn_bound))
    error('The number of bounds for measurement noise is not correct.');
end
if(length(state_bound)~=n||~isvector(state_bound))
    error('The number of bounds for states is not correct.');
end

% Check if the input and output are consistent.
if(length(input)~=length(output))
    error('The input length and output length are not consistent.');
end
% Check if the input is consistent with the model.
if(size(input,1)~=size(sys.mode.B,2))
    error('The input is not consistent with the model.');
end
% Check if the output is consistent with the model.
if(size(output,1)~=size(sys.mode.C,1))
    error('The output is not consistent with the model.');
end

% Construct the matrix M using system parameters.
% Construct the first part of M.
M1=[kron(eye(T),sys.mode.C) kron(eye(T),sys.Em) zeros(T*n_y,(T-1)*n)];
% Construct the second part of M.
M2_1=zeros((T-1)*n,T*n);
for i=1:T-1
    M2_1(n*(i-1)+1:n*i,n*(i-1)+1:n*(i+1))=[-sys.mode.A eye(n)];
end
M2_2=[zeros((T-1)*n,T*n_y) kron(eye(T-1),-sys.Ep)];
M2=[M2_1 M2_2];
% Combine M1 and M2 to obtain M.
M=[M1;M2];

% Construct the vector b using I/O data as well as the system parameters.
% Construct the first part of b.
b1=output(:,1:T)-sys.mode.D*input(:,1:T)-...
    bsxfun(@plus,zeros(n_y,T),sys.mode.f);
b1=reshape(b1,[],1);
% Construct the second part of b.
b2=sys.mode.B*input(:,1:T-1)+bsxfun(@plus,zeros(n,T-1),sys.mode.g);
b2=reshape(b2,[],1);
% Combine b1 and b2 to obtain b.
b=[b1;b2];

% Set up variables and constraints.
x=sdpvar(T*(2*n+n_y)-n,1);
constraints=[M*x==b];
for i=1:n
    if(state_bound(i)~=inf)
        constraints=[constraints,...
            norm(x(i:n:n*T),sys.state_norm(i))<=state_bound(i)];
    end
end
for i=1:n_y
    if(mn_bound(i)~=inf)
        constraints=[constraints,...
            norm(x(n*T+i:n_y:n*T+n_y*T),sys.mn_norm(i))<=mn_bound(i)];
    end
end
for i=1:n
    if(pn_bound(i)~=inf)
        constraints=[constraints,...
            norm(x((n+n_y)*T+i:n:end),sys.pn_norm(i))<=pn_bound(i)];
    end
end

% Check the feasibility.
options=sdpsettings('verbose',0,'solver',solver);
solution=optimize(constraints,[],options);
if(solution.problem==0)
    result=true;
else
    result=false;
end

end
