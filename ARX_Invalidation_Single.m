% This is a file running for system invalidation. The system model is an ARX
%   model with a single mode (i.e. only one submodel).

clear all

% Define system parameter.
A(:,:,1)=[1 2;2 2];
A(:,:,2)=[2 2;1 3];
C(:,:,1)=[2 1 2;3 1 5];
C(:,:,2)=[2 2 1;1 3 1];
f=[2;1];

% Set up noise bound.
pn_norm=[inf inf];
mn_norm=[inf inf];
pn_bound=[3 3];
mn_bound=[3 3];

% Creat the system model.
sys=ARXmodel(A,C,f,pn_norm,mn_norm);

degree=size(sys.mode.A,3);
n_in=size(C,2); % Input dimension.
n_y=size(A,1); % Output dimension.

% Set up desired input data and time horizon T (T may be greater than the length of input).
T=20;
input=[zeros(n_in,degree) ones(n_in,T-degree)];

% Set up initial conditions if needed.
ini_cond=[];

% Run simulation to obtain I/O data.
[y,~,~,~]=simulates(sys,input,T,ini_cond,pn_bound,mn_bound);

% Construct M and b of M*n=b for feasibility testing.
% Construct M.
M1=eye(n_y*(T-degree));
M2=zeros(n_y,n_y*T);
num_loop=T-degree;
buffer=eye(n_y);
for i=1:degree
    buffer(:,n_y*i+1)=sys.mode.A(:,:,i);
end
for i=1:num_loop
    M2(1+n_y*(i-1):n_y*i,1+n_y*(i-1):(degree+i)*n_y)=buffer;
end
M=[M1 M2];
% Construct b.
b=zeros(n_y,n_y*(T-degree));
for i=1:degree
    b=b+sys.mode.A(:,:,i)*y(:,degree+1-i:T-i)+sys.mode.C(:,:,i)*input(:,degree+1-i:T-i);
end
F=bsxfun(zeros(n_y,num_loop),f);
b=b-y(:,degree+1:T)+F;
b=reshape(b,[],1);

% Check if b is in the range of M with the constraint of bounded norm of n.
P=sdpavr((T-nA+1)*n_in+T*n_in,1);
C=[];
