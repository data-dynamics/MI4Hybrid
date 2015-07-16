function [y,p_noise,m_noise,switchseq]=swarx_sim(sys,input,ini_reg,...
    pn_bound,mn_bound,input_bound,switchseq,flag)

% This function simulates a switched ARX model with a giving input sequence.
%
% Arguments:
%   sys -- a user-defined system class (see StateSpace.m)
%   input -- input sequence (an n_i-by-T matrix) where n_i is the dimension
%            of input and T is the time horizon
%   ini_reg -- initial regressor for ARX models (an n_y-by-n matrix for n_y
%              outputs and order n)
%   pn_bound -- an n_y-D vector specifying the upper bound for process noise
%               where n is the number of outputs
%   mn_bound -- an n_y-D vector specifying the upper bound for measurement
%               noise where n_y is the number of outputs
%   input_bound -- an n_i-D vector specifying the upper bound for inputs
%                  where n_i is the number of inputs
%   switchseq -- a vector specifying the switching sequence that users want
%                to use
%   flag -- the simulation will begin with a fixed initial seed for the
%           random number generator if flag is set to 1
% Outputs:
%   y -- output of the system
%   p_noise -- the generated process noise
%   m_noise -- the generated measurement noise
%   switchseq -- the switching sequence
%
% Syntax:
%   [y,p_noise,m_noise,switchseq]=swarx_sim(sys,input);
%   [y,p_noise,m_noise,switchseq]=swarx_sim(sys,input,ini_cond);
%   [y,p_noise,m_noise,switchseq]=swarx_sim(sys,input,ini_cond,pn_bound,...
%                         mn_bound);
%   [y,p_noise,m_noise,switchseq]=swarx_sim(sys,input,ini_cond,pn_bound,...
%                         mn_bound,input_bound);
%   [y,p_noise,m_noise,switchseq]=swarx_sim(sys,input,ini_cond,pn_bound,...
%                         mn_bound,input_bound,switchseq);
%   [y,p_noise,m_noise,switchseq]=swarx_sim(sys,input,ini_cond,pn_bound,...
%                         mn_bound,input_bound,switchseq,flag);
%
% Author: Z. Luo, F. Harirchi and N. Ozay
% Date: July 15th, 2015

% Obtain model/input information.
num_arg=nargin;
n=size(sys.mode(1).A,3); % order of the system
n_y=size(sys.mode(1).A,1); % dimension of system output
n_i=size(sys.mode(1).C,2); % input dimension
N=size(sys.mode,2);   % number of modes

% Set up default values if parameters are not specified.
if(num_arg==2)
    ini_reg=zeros(n_y,n);
    pn_bound=zeros(n_y,1);
    mn_bound=zeros(n_y,1);
    input_bound=zeros(n_i,1)+inf;
    switchseq=[];
    flag=0;
elseif(num_arg==3)
    pn_bound=zeros(n_y,1);
    mn_bound=zeros(n_y,1);
    input_bound=zeros(n_i,1)+inf;
    switchseq=[];
    flag=0;
elseif(num_arg==5)
    input_bound=zeros(n_i,1)+inf;
    switchseq=[];
    flag=0;
elseif(num_arg==6)
    switchseq=[];
    flag=0;
elseif(num_arg==7)
    flag=0;
end

% Set up default values for empty arguments.
if(isempty(pn_bound))
    pn_bound=zeros(n_y,1);
end
if(isempty(mn_bound))
    mn_bound=zeros(n_y,1);
end
if(isempty(input_bound))
    input_bound=zeros(n_i,1)+inf;
end
if(isempty(flag))
    flag=0;
end

% Make the initial regressor be zero if it's not specified by users.
if(isempty(ini_reg))
    ini_reg=zeros(n_y,n);
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
if(length(input_bound)==1&&n_i>1)
    input_bound=ones(n_i,1)*input_bound;
    warning(['Input bound is a scalar, converted to a vector with '...
        'identical entries.']);
end

% Check the bounds.
if(length(pn_bound)~=n_y||~isvector(pn_bound))
    error('The number of bounds for process noise is not correct.');
end
if(length(mn_bound)~=n_y||~isvector(mn_bound))
    error('The number of bounds for measurement noise is not correct.');
end
if(length(input_bound)~=n_i||~isvector(input_bound))
    error('The number of bounds for inputs is not correct.');
end

% Check the input.
if(size(input,1)~=n_i)
    error('The input dimension is not correct.');
end
for i=1:n_i
    if(input_bound(i)~=inf)
        if(norm(input(i,:),sys.input_norm(i))>input_bound(i))
            warning(['At least one dimension of the input sequence '...
                'exceeds its bound']);
            break
        end
    end
end
T=size(input,2); % time horizon

% Use a random switching sequence if it is an empty input, otherwise use
% the user-specified switching sequence.
if(isempty(switchseq))
    switchseq=randi(N,[T-n,1]); % no need to have a length of n
end

% Creat process noise.
a=sys.pn_norm;
k=pn_bound;
p_noise=bounded_noise(a,k,T,flag);

% Creat measurement noise.
a=sys.mn_norm;
k=mn_bound;
m_noise=bounded_noise(a,k,T);

% Calculate the output using the switching sequence.
y=zeros(n_y,T); % pre-allocate memory
y(:,1:n)=ini_reg;
for i=n+1:T
    % Note: u(:,i) is not included when calculating y(:,i).
    % Note: p_noise(:,1:n) is not included for calculation.
    y(:,i)=zeros(n_y,1);
    for j=1:n
        y(:,i)=y(:,i)+sys.mode(switchseq(i-n)).A(:,:,j)*y(:,i-j)+...
            sys.mode(switchseq(i-n)).C(:,:,j)*input(:,i-j);
    end
    y(:,i)=y(:,i)+sys.mode(switchseq(i-n)).f+sys.Ep*p_noise(:,i);
end
y=y+sys.Em*m_noise;

end