function result=invalidation_arx(sys,input,output,pn_bound,mn_bound,solver)

% This function applies the invalidation algorithm for ARX models based on
% the system input and output. This invlidation algorithm is made to check
% the feasibility of a linear equation in the form of M*x=b where x is the
% variable vector constructed using noise vectors.
%
% Input/parameters:
%   sys -- the system model
%   input -- the input sequence
%   output -- the output sequence
%   pn_bound -- the bound for process noise
%   mn_bound -- the bound for measurement noise
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
%   result=invalidation_arx(sys,input,output,pn_bound,mn_bound);
%   result=invalidation_arx(sys,input,output,pn_bound,mn_bound,solver);
%
% Author: MI4Hybrid
% Date: June 22nd, 2015

% Check if the system model is valid for this function.
if(strcmp(sys.mark,'arx')~=1)
    error('The system model must be a non-switched ARX model.');
end

% Obtain the system model information.
T=size(input,2); % time horizon
n_y=size(sys.mode.A,1); % output dimension
degree=size(sys.mode.A,3); % degree of the system

% Use the default solver if it is not specified.
if(nargin==5)
    solver='cplex';
end

% Set up default values for empty paramters.
if(isempty(pn_bound))
    pn_bound=zeros(n_y,1);
end
if(isempty(mn_bound))
    mn_bound=zeros(n_y,1);
end
if(isempty(solver))
    solver='cplex';
end

% Check if the input and output are consistent.
if(length(input)~=length(output))
    error('The input length and output length are not consistent.');
end
% Check if the input is consistent with the model.
if(size(input,1)~=size(sys.mode.C,2))
    error('The input is not consistent with the model.');
end
% Check if the output is consistent with the model.
if(size(output,1)~=size(sys.mode.A,1))
    error('The output is not consistent with the model.');
end

% Convert scalars to vectors.
if(length(pn_bound)==1&&n_y>1)
    pn_bound=ones(n_y,1)*pn_bound;
    warning(['Bound for process noise is a scalar, converted to a '...
        'vector with identical entries.']);
end
if(length(mn_bound)==1&&n_y>1)
    mn_bound=ones(n_y,1)*mn_bound;
    warning(['Bound for measurement noise is a scalar, converted to'...
        ' a vector with identical entries.']);
end

% Check the bounds.
if(length(pn_bound)~=n_y||~isvector(pn_bound))
    error('The number of bounds for process noise is not correct.');
end
if(length(mn_bound)~=n_y||~isvector(mn_bound))
    error('The number of bounds for measurement noise is not correct.');
end

% Construct the matrix M using system parameters.
M1=kron(eye(T-degree),sys.Ep);
M2=zeros(n_y*(T-degree),n_y*T);
num_loop=T-degree;
buffer=zeros(n_y,n_y*(degree+1));
for i=1:degree
    buffer(:,n_y*(i-1)+1:n_y*i)=sys.mode.A(:,:,degree-i+1)*sys.Em;
end
buffer(:,n_y*degree+1:n_y*(degree+1))=-sys.Em;
for i=1:num_loop
    M2(1+n_y*(i-1):n_y*i,1+n_y*(i-1):(degree+i)*n_y)=buffer;
end
M=[-M1 M2];

% Construct the vector b using I/O data as well as the system parameters.
b=zeros(n_y,T-degree);
for i=1:degree
    b=b+sys.mode.A(:,:,i)*output(:,degree+1-i:T-i)+...
        sys.mode.C(:,:,i)*input(:,degree+1-i:T-i);
end
F=bsxfun(@plus,zeros(n_y,num_loop),sys.mode.f);
b=b-output(:,degree+1:T)+F;
b=reshape(b,[],1);

% Check the feasibility subject to the noise bound.
n=sdpvar((T-degree)*n_y+T*n_y,1);
constraints=[M*n==b];
for i=1:n_y
    if(pn_bound(i)~=inf)
        constraints=[constraints,...
            norm(n(i:n_y:n_y*(T-degree)),sys.pn_norm(i))<=pn_bound(i)];
    end
end
for i=1:n_y
    if(mn_bound(i)~=inf)
        constraints=[constraints,...
            norm(n(n_y*(T-degree)+i:n_y:end),sys.mn_norm(i))<=mn_bound(i)];
    end
end
options=sdpsettings('verbose',0,'solver',solver);
solution=optimize(constraints,[],options);
if(solution.problem==0)
    result=true;
else
    result=false;
end

end