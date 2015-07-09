function result=InvalidationARX(sys,input,output,pn_bound,mn_bound)

% This function applies the invalidation algorithm for ARX models based on
% the system input and output. This invlidation algorithm is made to check
% the feasibility of a linear equation in the form of M*x=b where x is the
% variable vector constructed using noise vectors.
%
% Input/parameters:
%   sys -- the system model
%   input -- the input sequence
%   output -- the output sequence
% Output:
%   result -- return "true" if the system is validated, otherwise "false"
%
% Author: MI4Hybrid
% Date: June 22nd, 2015

% Check if the system model is valid for this function.
if(strcmp(sys.mark,'arx')~=1)
    error('The system model must be a non-switched ARX model.');
end
% Check if the input and output are consistent.
if(length(input)~=length(output))
    error('The input length and output length are not consistent.');
end
% Check if the input is consistent with with the model.
if(size(input,1)~=size(sys.mode.C,2))
    error('The input is not consistent with the model.');
end
% Check if the output is consistent with the model.
if(size(output,1)~=size(sys.mode.A,1))
    error('The output is not consistent with the model.');
end

% Obtain the system model information.
T=size(input,2); % time horizon
n_y=size(sys.mode.A,1); % output dimension
degree=size(sys.mode.A,3); % degree of the system

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
F=bsxfun(@plus,zeros(n_y,num_loop),sys.f);
b=b-output(:,degree+1:T)+F;
b=reshape(b,[],1);

% Check the feasibility subject to the noise bound.
n=sdpvar((T-degree)*n_y+T*n_y,1);
constraints=[M*n==b];
for i=1:n_y
    constraints=[constraints,...
        norm(n(i:n_y:n_y*(T-degree)),sys.pn_norm(i))<=pn_bound(i)];
end
for i=1:n_y
    constraints=[constraints,...
        norm(n(n_y*(T-degree)+i:n_y:end),sys.mn_norm(i))<=mn_bound(i)];
end
options=sdpsettings('verbose',0,'solver','mosek');
solution=optimize(constraints,[],options);
if(solution.problem==0)
    result=true;
else
    result=false;
end

end