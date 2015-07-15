function [y,p_noise,m_noise,switchseq]=simulates(sys,input,ini_cond,...
    pn_bound,mn_bound,input_bound,state_bound,switchseq,flag)

% This function can simulate (switched) state-space models, (switched) ARX
% models, and (not switched) polynomial models.
%
% Inputs/Parameters:
%   sys -- a user-defined system class (see StateSpace.m, ARXmodel.m, and
%          polymodel.m)
%   input -- input sequence (an n_i-by-T matrix) where n_i is the dimension
%            of input and T is its length
%   ini_cond -- initial conditions for state-space models (an n-by-1 column
%               vector for n states) or initial regressors for ARX models (
%               an n-by-1 column vector for order n)
%   pn_bound -- an n_p-dimensional vector specifying the bound for process
%               noise where n_p is the dimension of process noise
%   mn_bound -- an n_m-D vector specifying the bound for process noise
%               where n_m is the dimension of measurement noise
%   switchseq -- switching sequence for switched models that can be
%                specified by users (empty if not specified)
%   flag -- the simulation will begin with a fixed initial seed for the
%           random number generator if flag is set to 1
% Output:
%   y -- output of the system
%   p_noise -- the generated process noise
%   m_noise -- the generated measurement noise
%   switchseq -- the switching sequence used for simulation
%
% Syntax:
%   [y,p_noise,m_noise,switchseq]=simulates(sys,input);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,input,ini_cond);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,input,ini_cond,pn_bound,...
%                         mn_bound);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,input,ini_cond,pn_bound,...
%                         mn_bound,input_bound);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,input,ini_cond,pn_bound,...
%                         mn_bound,input_bound,state_bound);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,input,ini_cond,pn_bound,...
%                         mn_bound,input_bound,state_bound,switchseq);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,input,ini_cond,pn_bound,...
%                         mn_bound,input_bound,state_bound,switchseq,flag);
%
% Author: Z. Luo, F. Harirchi and N. Ozay
% Date: May 22nd, 2015

% Set up default values if parameters are not specified.
if(nargin==2)
    ini_cond=[];
    pn_bound=zeros(size(sys.mode(1).A,1),1);
    mn_bound=zeros(size(sys.mode(1).C,1),1);
    input_bound=
    state_bound=
    flag=0;
    switchseq=[];
end
if(nargin==4)
    pn_bound=zeros(size(sys.mode(1).A,1),1);
    mn_bound=zeros(size(sys.mode(1).C,1),1);
    flag=0;
    switchseq=[];
end
if(nargin==6)
    flag=0;
    switchseq=[];
end
if(nargin==7)
    switchseq=[];
end

% Convert a scalar (pn_bound or mn_bound) to a vector with identical entries.
if(length(pn_bound)==1&&size(sys.mode(1).A,1)>1)
    pn_bound=zeros(size(sys.mode(1).A,1),1)+pn_bound;
    warning('Input bound of process noise is a scalar, converted to a vector with identical entries.');
end
if(length(mn_bound)==1&&size(sys.mode(1).C,1)>1)
    mn_bound=zeros(size(sys.mode(1).C,1),1)+mn_bound;
    warning('Input bound of measurement noise is a scalar, converted to a vector with identical entries.');
end

switch sys.mark
    case 'swss'
        [y,p_noise,m_noise,switchseq]=swss_sim(sys,ini_cond,input,T,switchseq,pn_bound,mn_bound,flag);
    case 'ss'
        [y,p_noise,m_noise,switchseq]=ss_sim(sys,ini_cond,input,T,pn_bound,mn_bound,flag);
    case 'swarx'
        [y,p_noise,m_noise,switchseq]=swarx_sim(sys,ini_cond,input,T,switchseq,pn_bound,mn_bound,flag);
    case 'arx'
        [y,p_noise,m_noise,switchseq]=arx_sim(sys,ini_cond,input,T,pn_bound,mn_bound,flag);
    case 'poly'
        [y,p_noise,m_noise,switchseq]=poly_sim(sys,ini_cond,input,T,pn_bound,mn_bound,flag);
end

end

%--------------Simulate Switched State-Space Model--------------
function [y,p_noise,m_noise,switchseq]=swss_sim(sys,ini_cond,input,T,switchseq,pn_bound,mn_bound,flag)

% Use a random switching sequence if it is an empty argument, otherwise use
%   the user-specified switching sequence.
if(isempty(switchseq))
    N=size(sys.mode,2);   % number of modes
    switchseq=randi(N,[T,1]);
end

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

% Calculate the output using the switching sequence.
n=size(sys.mode(1).C,1); % dimension of system output
y=zeros(n,T); % pre-allocate memory
x=ini_cond;
for i=1:T
    x=sys.mode(switchseq(i)).A*x+sys.mode(switchseq(i)).B*input(:,i)+sys.g(:,switchseq(i))+sys.Ep*p_noise(:,i);
    y(:,i)=sys.mode(switchseq(i)).C*x+sys.mode(switchseq(i)).D*input(:,i)+sys.f(:,switchseq(i));
end
y=y+sys.Em*m_noise;

end

%--------------Simulate Switched ARX Model--------------
function [y,p_noise,m_noise,switchseq]=swarx_sim(sys,ini_cond,input,T,switchseq,pn_bound,mn_bound,flag)

n=size(sys.mode(1).A,3); % order of the system
out_dim=size(sys.mode(1).A,1); % dimension of system output

% Use a random switching sequence if it is an empty input, otherwise use
%   the user-specified switching sequence.
if(isempty(switchseq))
    N=size(sys.mode,2);   % number of modes
    switchseq=randi(N,[T-n,1]);
end

% Make the initial condition be zero if it's not specified by users.
if(isempty(ini_cond))
    ini_cond=zeros(out_dim,n);
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
y=zeros(out_dim,T); % pre-allocate memory
y(:,1:n)=ini_cond;
for i=n+1:T
    % Note: u(:,i) is not included when calculating y(:,i).
    y(:,i)=zeros(out_dim,1);
    for j=1:n
        y(:,i)=y(:,i)+sys.mode(switchseq(i-n)).A(:,:,j)*y(:,i-j)+sys.mode(switchseq(i-n)).C(:,:,j)*input(:,i-j);
    end
    y(:,i)=y(:,i)+sys.f(:,switchseq(i))+sys.Ep*p_noise(:,i); % add process noise
end
y=y+sys.Em*m_noise; % add measurement noise

end

%--------------Simulate State-Space Model--------------
function [y,p_noise,m_noise,switchseq]=ss_sim(sys,ini_cond,input,T,pn_bound,mn_bound,flag)

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
n=size(sys.mode.C,1); % dimension of output
y=zeros(n,T); % pre-allocate memory
x=ini_cond;
for i=1:T
    y(:,i)=sys.mode.C*x+sys.mode.D*input(:,i)+sys.f;
    x=sys.mode.A*x+sys.mode.B*input(:,i)+sys.g+sys.Ep*p_noise(:,i);
end
y=y+sys.Em*m_noise;

% The (trivial) switching sequence.
switchseq=ones(T,1);

end

%--------------Simulate ARX Model--------------
function [y,p_noise,m_noise,switchseq]=arx_sim(sys,ini_cond,input,T,pn_bound,mn_bound,flag)

out_dim=size(sys.mode.A,1); % dimension of system output
n=size(sys.mode.A,3); % order of the system

% Make the initial condition be zero if it's not specified by users.
if(isempty(ini_cond))
    ini_cond=zeros(out_dim,n);
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
y=zeros(out_dim,T); % pre-allocate memory
y(:,1:n)=ini_cond;
for i=n+1:T
    % Note: u(:,i) is not included when calculating y(:,i).
    y(:,i)=zeros(out_dim,1);
    for j=1:n
        y(:,i)=y(:,i)+sys.mode.A(:,:,j)*y(:,i-j)+sys.mode.C(:,:,j)*input(:,i-j);
    end
    y(:,i)=y(:,i)+sys.f+sys.Ep*p_noise(:,i); % add process noise
end
y=y+sys.Em*m_noise; % add measurement noise

% The (trivial) switching sequence.
switchseq=ones(T-n,1);

end

%--------------Simulate Polynomial Model--------------
function [y,p_noise,m_noise,switchseq]=poly_sim(sys,ini_cond,input,T,pn_bound,mn_bound,flag)

n=size(sys.coeffmat,1); % number of states/outputs
n_i=size(sys.degmat,2)-n; % input dimension

% Make the initial condition be zero if it's not specified by users.
if(isempty(ini_cond))
    ini_cond=zeros(out_dim,n);
end

% not finished yet






end