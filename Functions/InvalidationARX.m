function result=InvalidationARX(sys,input,output)

% This function applies the invalidation algorithm for ARX models based on
% the system input and output.
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

switch sys.mark
    case 'arx'
        result=InvalidationNARX(sys,input,output);
    case 'swarx'
        result=InvalidationSARX(sys,input,output);
end

end

%--------------Invalidate Non-switched ARX Model--------------
function result=InvalidationNARX(sys,input,output)

% This invlidation algorithm checks the feasibility of a linear equation in
% the form of M*x=b where x is the variable vector.

% It is possible that the input and output have different lengths since the
% function "simulates" will pad zeros to the input if the length of input
% is less than the specified time horizon T.
if(length(input)<length(output))
    %pad zero?
elseif(length(input)>length(output))
    %cut off?
end

% Construct the matrix M.
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

% Construc the vector b.
for i=1:degree
    b=b+sys.mode.A(:,:,i)*A_factor*y(:,degree+1-i:T-i)+sys.mode.C(:,:,i)*C_factor*input(:,degree+1-i:T-i);
end
F=bsxfun(@plus,zeros(n_y,num_loop),f);
b=b-y(:,degree+1:T)+F;
b=reshape(b,[],1);


% Check the feasibility subject to the noise bound.


end

%--------------Invalidate Switched ARX Model--------------
function result=InvalidationSARX(sys,input,output)

% Dummy implementation
result=true;

end