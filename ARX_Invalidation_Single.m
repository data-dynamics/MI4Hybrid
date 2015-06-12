% This is a file running for system invalidation. The system model is an ARX
%   model with a single mode (i.e. only one submodel).

% Define system parameter.
A=[];
C=[];
f=[];

% Set up noise bound.
pn_norm=[];
mn_norm=[];
pn_bound=[];
mn_bound=[];

% Creat the system model.
sys=ARXmodel();

% Set up desired input data and time horizon T (T may be greater than the length of input).
input=[];
T=0;

% Set up initial conditions if needed.
ini_cond=[];

% Run simulation to obtain I/O data.
[y,~,~,~]=simulates(sys,ini_cond,input,T,pn_bound,mn_bound);

% Apply system invalidation.
