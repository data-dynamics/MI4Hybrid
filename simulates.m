function [y,p_noise,m_noise,switchseq]=simulates(sys,ini_cond,input,T,pn_bound,mn_bound,flag,switchseq)

% This function can simulate (switched) state-space models, (switched) ARX models, and (switched)
%   polynomial models.
%
% Inputs/Parameters:
%   sys -- a user-defined system class (see StateSpace.m and ARXmodel.m)
%   ini_cond -- initial conditions for state-space models (an n-by-1 column vector for n states)
%               or initial regressors for ARX models (an n-by-1 column vector for order n)
%   input -- input sequence (an T-by-1 column vector)
%   T -- time horizon
%   pn_bound -- an n_p-D vector specifying the bound for process noise where n_p is the dimension
%               of process noise
%   mn_bound -- an n_m-D vector specifying the bound for process noise where n_m is the dimension
%               of measurement noise
%   switchseq -- switching sequence (for switched models) that can be specified (empty if not
%                specified)
%   flag -- the simulation will begin with a fixed initial seed for the random number generator
%           if flag is set to 1
% Output:
%   y -- output of the system
%   p_noise -- the generated process noise
%   m_noise -- the generated measurement noise
%   switchseq -- the switching sequence used for simulation
%
% Syntax:
%   [y,p_noise,m_noise,switchseq]=simulates(sys,ini_cond,input,T);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,ini_cond,input,T,pn_bound,mn_bound);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,ini_cond,input,T,pn_bound,mn_bound,flag);
%   [y,p_noise,m_noise,switchseq]=simulates(sys,ini_cond,input,T,pn_bound,mn_bound,flag,switchseq);
%
% Author: MI4Hybrid
% Date of Last Modification: May 22nd, 2015

% Set up default values if parameters are not specified.
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

% Add zeros to the end of input if T is greater than the length of input.
[input_dim,input_len]=size(input);
if(T>input_len)
    input=[input zeros(input_dim,T-input_len)];
    warning('Input length is less than time horizon, added zeros to the end of input.');
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
end

end

%--------------Simulate Switched State-Space Model--------------
function [y,p_noise,m_noise,switchseq]=swss_sim(sys,ini_cond,input,T,switchseq,pn_bound,mn_bound,flag)

% Use a random switching sequence if it is an empty input, otherwise use
%   the user-specified switching sequence.
if(isempty(switchseq))
    N=size(sys.mode,2);   % number of modes
    switchseq=randi(N,[T,1]);
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
y=zeros(n,T);
x=ini_cond;
for i=1:T
    x=sys.mode(switchseq(i)).A*x+sys.mode(switchseq(i)).B*input(:,i)+sys.g(:,switchseq(i))+sys.Ep*p_noise(:,i);
    y(:,i)=sys.mode(switchseq(i)).C*x+sys.mode(switchseq(i)).D*input(:,i)+sys.f(:,switchseq(i));
end
y=y+sys.Em*m_noise;

end

%--------------Simulate Switched ARX Model--------------
function [y,p_noise,m_noise,switchseq]=swarx_sim(sys,ini_cond,input,T,switchseq,pn_bound,mn_bound,flag)

% Use a random switching sequence if it is an empty input, otherwise use
%   the user-specified switching sequence.
if(isempty(switchseq))
    N=size(sys.mode,2);   % number of modes
    n=size(sys.mode(1).A,2); % order of the system
    switchseq=randi(N,[T-n,1]);
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
out_dim=size(sys.mode(1).A,1); % dimension of system output
y=zeros(out_dim,T);
m=size(sys.mode(1).C,2);
y(:,1:n)=ini_cond;
% If m>n, return error or pad zero?
for i=n+1:T
    % Note: u(:,i) is not included when calculating y(:,i).
    y(:,i)=sys.mode(switchseq(i-n)).A*y(:,i-n:i-1)'+sys.mode(switchseq(i-n)).C*input(:,i-m:i-1)';
end
y(:,i)=y(:,i)+sys.Ep*p_noise; % add process noise
y(:,i)=y(:,i)+sys.Em*m_noise; % add measurement noise

end

%--------------Simulate State-Space Model--------------
function [y,p_noise,m_noise,switchseq]=ss_sim(sys,ini_cond,input,T,pn_bound,mn_bound,flag)

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
y=zeros(n,T);
x=ini_cond;
for i=1:T
    x=sys.mode.A*x+sys.mode.B*input(:,i)+sys.Ep*p_noise(:,i);
    y(:,i)=sys.mode.C*x+sys.mode.D*input(:,i);
end
y=y+sys.Em*m_noise;

% The (trivial) switching sequence.
switchseq=ones(T,1);

end

%--------------Simulate ARX Model--------------
function [y,p_noise,m_noise,switchseq]=arx_sim(sys,ini_cond,input,T,pn_bound,mn_bound,flag)

% Creat process noise.
a=sys.pn_norm;
k=pn_bound;
p_noise=bounded_noise(a,k,T,flag);

% Creat measurement noise.
a=sys.mn_norm;
k=mn_bound;
m_noise=bounded_noise(a,k,T);

% Calculate the output.
out_dim=size(sys.mode.A,1); % dimension of system output
y=zeros(out_dim,T);
n=size(sys.mode.A,2); % order of the system
m=size(sys.mode.C,2);
y(:,1:n)=ini_cond;
% If m>n, return error or pad zero?
for i=n+1:T
    % Note: u(:,i) is not included when calculating y(:,i).
    y(:,i)=sys.mode.A*y(:,i-n:i-1)'+sys.mode.C*input(:,i-m:i-1)';
end
y(:,i)=y(:,i)+sys.Ep*p_noise; % add process noise
y(:,i)=y(:,i)+sys.Em*m_noise; % add measurement noise

% The (trivial) switching sequence.
switchseq=ones(T,1);

end