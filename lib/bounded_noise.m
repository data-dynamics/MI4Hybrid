function noise=bounded_noise(a,k,T,flag)

% This function creats vector valued (scalar valued if the dimension of the
% noise is 1) noise whose norm is bounded by k. Each dimension of the noise
% is uniformly distributed.
%
% Inputs:
%   a -- an n-by-1 vector corresponding to the norm types (for n dimension)
%   k -- an n-by-1 vector corresponding to the upper bounds (for n dimension)
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
% Date: May 22nd, 2015

% Change a scalar bound to a vector having the same entries.
if(length(k)==1&&length(a)~=1)
    k=zeros(length(a),1)+k;
    warning(['Input bound is a scalar, converted to a vector having the'...
        ' same entries.']);
end

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
        noise(i,:)=(rand(1,T)-0.5)*2;
        noise_norm=norm(noise(i,:),a(i));
        noise(i,:)=noise(i,:)*k(i)/noise_norm;
    else
        noise(i,:)=(rand(1,T)-0.5)*2*k(i);
    end
end

end