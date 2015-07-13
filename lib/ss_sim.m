function [y,p_noise,m_noise]=ss_sim(sys,input,ini_cond,pn_bound,...
    mn_bound,input_bound,state_bound,flag)

% This function simulates a non-switched state-space model with a giving
% input sequence.
%
% Arguments:
%   sys -- a user-defined system class (see StateSpace.m)
%   input -- input sequence (an n_i-by-T matrix) where n_i is the dimension
%            of input and T is the time horizon
%   ini_cond -- initial condition for state-space models (an n-by-1 column
%               vector for n states)
%   pn_bound -- an n-D vector specifying the upper bound for process noise
%               where n is the number of states
%   mn_bound -- an n_y-D vector specifying the upper bound for measurement
%               noise where n_y is the number of outputs
%   input_bound -- an n_i-D vector specifying the upper bound for inputs
%                  where n_i is the number of inputs
%   state_bound -- an n-D vector specifying the upper bound for states
%                  where n is the number of states
%   flag -- the simulation will begin with a fixed initial seed for the
%           random number generator if flag is set to 1
% Outputs:
%   y -- output of the system
%   p_noise -- the generated process noise
%   m_noise -- the generated measurement noise
%
% Syntax:
%   [y,p_noise,m_noise]=ss_sim(sys,input);
%   [y,p_noise,m_noise]=ss_sim(sys,input,ini_cond);
%   [y,p_noise,m_noise]=ss_sim(sys,input,ini_cond,pn_bound,mn_bound);
%   [y,p_noise,m_noise]=ss_sim(sys,input,ini_cond,pn_bound,mn_bound,...
%                             input_bound);
%   [y,p_noise,m_noise]=ss_sim(sys,input,ini_cond,pn_bound,mn_bound,...
%                             input_bound,state_bound);
%   [y,p_noise,m_noise]=ss_sim(sys,input,ini_cond,pn_bound,mn_bound,...
%                             input_bound,state_bound,flag);
%
% Author: Z. Luo, F. Harirchi and N. Ozay
% Date: July 13th, 2015

% Set up default values.
num_arg=nargin;
n=size(sys.mode.A,1); % dimension of state
n_y=size(sys.mode.C,1); % dimension of output
n_i=size(sys.mode.B,2); % dimension of input
if(num_arg==2)
    ini_cond=[];
    pn_bound=zeros(n,1);
    mn_bound=zeros(n_y,1);
    input_bound=zeros(n_i,1)+inf;
    state_bound=zeros(n,1)+inf;
    flag=0;
elseif(num_arg==3)
    pn_bound=zeros(n,1);
    mn_bound=zeros(n_y,1);
    input_bound=zeros(n_i,1)+inf;
    state_bound=zeros(n,1)+inf;
    flag=0;
elseif(num_arg==5)
    input_bound=zeros(n_i,1)+inf;
    state_bound=zeros(n,1)+inf;
    flag=0;
elseif(num_arg==6)
    state_bound=zeros(n,1)+inf;
    flag=0;
elseif(num_arg==7)
    flag=0;
end

% Check the input.
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

% Make the initial condition be zero if it's not specified by users.
if(isempty(ini_cond))
    ini_cond=zeros(size(sys.mode(1).A,1),1);
end

% Creat process noise.
a=sys.pn_norm;
k=pn_bound;
p_noise=bounded_noise(a,k,T,flag);

% Creat measurement noise.
a=sys.mn_norm;
k=mn_bound;
m_noise=bounded_noise(a,k,T);

% Calculate the output.
y=zeros(n_y,T); % pre-allocate memory
x(:,1)=ini_cond;
for i=1:T
    y(:,i)=sys.mode.C*x(:,i)+sys.mode.D*input(:,i)+sys.f;
    x(:,i+1)=sys.mode.A*x(:,i)+sys.mode.B*input(:,i)+sys.g...
        +sys.Ep*p_noise(:,i);
end
y=y+sys.Em*m_noise;

% Check the state.
for i=1:n
    if(state_bound(i)~=inf)
        if(norm(x(i,1:T),sys.state_norm(i))>state_bound(i))
            warning(['At least one dimension of the state sequence '...
                'exceeds its bound']);
            break
        end
    end
end

end