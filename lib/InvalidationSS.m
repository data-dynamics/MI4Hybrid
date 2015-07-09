function result=InvalidationSS(sys,input,output,pn_bound,mn_bound)

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
% Output:
%   result -- return "true" if the system is validated, otherwise "false"
%
% Syntax:
%   result=InvalidationSS(sys,input,output);
%   result=InvalidationSS(sys,input,output,pn_bound,mn_bound);
%
% Author: MI4Hybrid
% Date: July 8th, 2015

% Use the default bounds if they are not specified.
if(nargin==3&&~isnan(sys.pn_bound)&&~isnan(sys.mn_bound))
    pn_bound=sys.pn_bound;
    mn_bound=sys.mn_bound;
elseif(nargin~=5)
    error('Noise bounds are not specified.');
end

% Check if the system model is valid for this function.
if(strcmp(sys.mark,'ss')~=1)
    error('The system model must be a non-switched state-space model.');
end
% Check if the input and output are consistent.
if(length(input)~=length(output))
    error('The input length and output length are not consistent.');
end
% Check if the input is consistent with with the model.
if(size(input,1)~=size(sys.mode.B,2))
    error('The input is not consistent with the model.');
end
% Check if the output is consistent with the model.
if(size(output,1)~=size(sys.mode.C,1))
    error('The output is not consistent with the model.');
end

% Obtain system model information.
T=size(input,2); % time horizon
n=size(sys.mode.C,2); % state dimension
n_y=size(sys.mode.C,1); % output dimension

% Construct the matrix M using system parameters.
% Construct the first part of M.
buffer=sys.mode.C*sys.mode.A;
M1=kron(eye(T-1),buffer);
buffer=zeros((T-1)*n_y,n);
M1=[M1 buffer];
buffer=zeros(n_y,T*n);
buffer(:,1:n)=sys.mode.C;
M1=[buffer;M1];
buffer=zeros(n_y,T*n);
buffer(:,n*(T-1)+1:n*T)=sys.mode.C;
M1=[M1;buffer];
% Construct the second part of M.
M2=zeros(n_y*(T+1),n*T);
buffer=sys.mode.C*sys.Ep;
M2(n_y+1:n_y*T,1:n*(T-1))=kron(eye(T-1),buffer);
% Construct the third part of M.
M3=kron(eye(T),sys.Em);
buffer=zeros(n_y,n_y*T);
buffer(:,n_y*(T-1)+1:n_y*T)=sys.Em;
M3=[M3;buffer];
% Combine M1, M2, M3 to obtain M.
M=[M1 M2 M3];

% Construct the vector b using I/O data as well as the system parameters.
b=output(:,2:end)-...
    sys.mode.C*sys.mode.B*input(:,1:T-1)-...
    sys.mode.D*input(:,2:end)-...
    sys.mode.C*bsxfun(@plus,zeros(n,T-1),sys.g)-...
    bsxfun(@plus,zeros(n_y,T-1),sys.f);
b=reshape(b,[],1);
b=[output(:,1)-sys.mode.D*input(:,1)-sys.f;b];
b=[b;output(:,T)-sys.mode.D*input(:,T)-sys.f];

% Set up variables and constraints.
x=sdpvar(T*(2*n+n_y),1);
constraints=[M*x==b];
for i=1:n
    constraints=[constraints,...
        norm(x(i:n:n*T),sys.state_norm(i))<=sys.state_bound(i)];
end
for i=1:n
    constraints=[constraints,...
        norm(x(n*T+i:n:n*2*T),sys.pn_norm(i))<=pn_bound(i)];
end
for i=1:n_y
    constraints=[constraints,...
        norm(x(2*n*T+i:n_y:end),sys.mn_norm(i))<=mn_bound(i)];
end

% Check the fesibility.
options=sdpsettings('verbose',0,'solver','mosek');
solution=optimize(constraints,[],options);
if(solution.problem==0)
    result=true;
else
    result=false;
end

end