function noise=bounded_noise(a,k,T,flag)

% This function creats vector valued (scalar valued if the dimension of the
%   noise vectors is 1) noise whose a1-a2 norm is bounded by k.
%
% Inputs:
%   a -- an n-by-1 vector corresponding to the norm types
%   k -- an n-by-1 vector corresponding to the upper bounds
%   T -- time horizon
%   flag -- set to 1 if the initial seed needs to be reset to a fixed value
% Outputs:
%   noise -- the output noise
%
% Syntax:
%   noise=bounded_noise(a,k,T);
%   noise=bounded_noise(a,k,T,flag);
%
% Author: MI4Hybrid
% Date of Last Modification: May 22nd, 2015

% Find the dimension of noise.
if(length(a)~=length(k))
    error('The number of norm types must match with the number of bounds.');
else
    n=length(a); % noise dimension
end

% Check if the norm types are valid.
if(~isempty(find(a<1,1)))
    error('For an l_p norm, p must be greater than or equal to 1.');
end

% The value of flag is 0 by default.
if(nargin==3)
    flag=0;
end

% Reset the initial seed to a fixed default value if flag is 1.
if(flag==1)
    rng('default');
end

% Creat noise.
noise=zeros(n,T);
for i=1:n
    if(a(i)>=1&&a(i)~=inf)
        element_bound=nthroot(k(i)^a(i)/T,a(i));
        noise(i,:)=(rand(1,T)-0.5)*2*element_bound;
    else
        noise(i,:)=(rand(1,T)-0.5)*2*k(i);
    end
end

end