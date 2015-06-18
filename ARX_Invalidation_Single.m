% This is a file running for system invalidation. The system model is an ARX
%   model with a single mode (i.e. only one submodel).

clear all

A_factor=1;
C_factor=1;

% Define system parameter.
A(:,:,1)=[0.4747 0.0628;-0.3424 1.2250];
A(:,:,2)=[0.5230 -0.1144;0.3574 -0.2513];
C(:,:,1)=[1 0 0;0 0 0];
C(:,:,2)=[1 0 0;0 0 0];
f=[0;0];

% Set up noise bound.
pn_norm=[inf inf];
mn_norm=[inf inf];
pn_bound=[2 2];
mn_bound=[2 2];

% Creat the system model.
sys=ARXmodel(A,C,f,pn_norm,mn_norm);

degree=size(sys.mode.A,3);
n_in=size(sys.mode.C,2); % Input dimension.
n_y=size(sys.mode.A,1); % Output dimension.

% Set up desired input data and time horizon T (T may be greater than the length of input).
T=500;
input=[zeros(n_in,degree) ones(n_in,T-degree)];

% Set up initial conditions if needed.
ini_cond=[];

% Run simulation to obtain I/O data.
[y,p_noise,m_noise,~]=simulates(sys,input,T,ini_cond,pn_bound,mn_bound);

% Construct M and b of M*n=b for feasibility testing.
% Construct M.
M1=eye(n_y*(T-degree));
M2=zeros(n_y*(T-degree),n_y*T);
num_loop=T-degree;
buffer=zeros(n_y,n_y*(degree+1));
for i=1:degree
    buffer(:,n_y*(i-1)+1:n_y*i)=sys.mode.A(:,:,degree-i+1)*A_factor;
end
buffer(:,n_y*degree+1:n_y*(degree+1))=-eye(n_y);
for i=1:num_loop
    M2(1+n_y*(i-1):n_y*i,1+n_y*(i-1):(degree+i)*n_y)=buffer;
end
M=[-M1 M2];
% Construct b.
b=zeros(n_y,T-degree);
for i=1:degree
    b=b+sys.mode.A(:,:,i)*A_factor*y(:,degree+1-i:T-i)+sys.mode.C(:,:,i)*C_factor*input(:,degree+1-i:T-i);
end
F=bsxfun(@plus,zeros(n_y,num_loop),f);
b=b-y(:,degree+1:T)+F;
b=reshape(b,[],1);

% Check if b is in the range of M with the constraint of bounded norm of n.
n=sdpvar((T-degree)*n_y+T*n_y,1);
Cons=[M*n==b, norm(n,inf)<=2];
tic
options =sdpsettings('verbose',1,'solver','mosek');
solution=optimize(Cons,[],options);
toc
if(solution.problem==0)
    disp('Valid');
elseif(solution.problem==1)
    disp('Invalid');
end