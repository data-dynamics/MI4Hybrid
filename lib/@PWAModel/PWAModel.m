classdef PWAModel
    % Beta version: no process noise and measurement noise considered
    %
    % The class represents a discrete-time piecewise swith affine(PWA) 
    % system model. A general PWA system with a switching
    % sequence sigma has the following form:
    %   x[k+1] = A[sigma[k]]*x[k] + B[sigma[k]]*u[k] + f[sigma[k]] + Ep*pn[k]
    %   y[k] = C[sigma[k]]*x[k] + D[sigma[k]]*u[k] + g[sigma[k]]+ Em*mn[k]
    %   sigma[k] = { i | P(i)*x[k] + M(i) <=0 }
    % where pn is the process noise and mn is the measurement noise.
    % Notice that the last nonpositive condition means each element of 
    % vector P(i)*x[k] + M(i) is nonpositive. 
    % 
    % Constructor syntax:
    %   sys=PWAModel(A,B,C,D,P,M);
    %   sys=PWAModel(A,B,C,D,P,M,f);
    %   sys=PWAModel(A,B,C,D,P,M,f,g);
    %   sys=PWAModel(A,B,C,D,P,M,f,g,pn_norm, mn_norm);
    %   sys=PWAModel(A,B,C,D,P,M,f,g,pn_norm, mn_norm, Ep, Em);
    %   sys=PWAModel(A,B,C,D,P,M,f,g,pn_norm, mn_norm, Ep, Em,input_norm);
    %   sys=PWAModel(A,B,C,D,P,M,f,g,pn_norm, mn_norm, Ep, Em,input_norm,state_norm);
    %
    % Author: Z. Luo, F. Harirchi and N. Ozay
    % Modified: J. Liu
    % Date: Sep 14th, 2016

    % Notations:
    %   n -- number of states
    %   n_y -- number of outputs
    %   n_i -- number of inputs
    %   n_mode -- number of modes
    %   n_hyper -- number of hyperplanes to consider activated mode
    
    properties(SetAccess = protected)
        % A set of discrete-time state-space modes.
        % e.g. mode(i).A, mode(i).B, mode(i).C, mode(i).D, mode(i).f, 
        % mode(i).g, mode(i).P and mode(i).M represent the i-th mode.
        mode
        % An n-by-1 column vector representing the norm types of process
        % noise.
        pn_norm
        % An n_y-by-1 column vector representing the norm types of
        % measurement noise.
        mn_norm
        % Ep is an n-by-n matrix.
        Ep
        % Em is an n_y-by-n_y matrix.
        Em
        % An n-by-1 column vector representing the norm types of states.
        state_norm
        % An n_i-by-n_i column vector representing the norm types of inputs.
        input_norm
        % A mark (a string) stating that the model is a PWA system.
        mark
    end
    
    methods
        
        % If there is only one mode, the model is not switchable.
        function sys=PWAModel(A,B,C,D,P,M,f,g,pn_norm, mn_norm, Ep, Em,input_norm,state_norm)
            
            % Check A, B, C, D, f, g, P, M.
            if(size(A,3)~=size(B,3)||size(A,3)~=size(C,3)||size(A,3)~=size(D,3)...
                    ||size(A,3)~=size(P,3))
                error('A, B, C, D, P must have the same number of matrices.');
            elseif(size(A,1)~=size(A,2)||size(A,1)~=size(B,1))
                error(['A must be a (or a set of) square matrices, or A '...
                      'and B not consistent.']);
            elseif(size(C,2)~=size(A,1)||size(C,1)~=size(D,1))
                error('C is not consistent with A or D.');
            elseif(size(P,1)~=size(M,1))
                error('P and M not consistent');
            elseif(size(P,2)~=size(A,1))
                error('P and A not consistent');
            elseif(size(M,2)~=size(A,3))
                error('M must have the same number of rows as the number of modes');
            else
                n_mode=size(A,3); % number of modes
                n_y=size(C,1); % number of outputs
                n=size(A,1); % number of states
                n_i=size(B,2); % number of inputs
                n_hyper = size(P,1); % number of hyperplanes
            end
            
            
            
            % Indicate the model type.
            if(n_mode>=1)
                sys.mark='pwa';
            else
                error('There should be at least one mode.');
            end
            
            % Set up default values if parameters are not specified.
            if(nargin==6)
                f=zeros(n,n_mode);
                g=zeros(n_y,n_mode);
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n);
                Em=eye(n_y);
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==7)
                g=zeros(n_y,n_mode);
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n);
                Em=eye(n_y);
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==8)
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n);
                Em=eye(n_y);
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==10)
                Ep=eye(n);
                Em=eye(n_y);
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==12)
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==13)
                state_norm=zeros(n,1)+inf;
            end
            
            % Set up default values for empty inputs.
            if(isempty(g))
                g=zeros(n_y,n_mode);
            end
            if(isempty(f))
                f=zeros(n,n_mode);
            end
            if(isempty(pn_norm))
                pn_norm=zeros(n,1)+inf;
            end
            if(isempty(mn_norm))
                mn_norm=zeros(n_y,1)+inf;
            end
            if(isempty(Ep))
                Ep=eye(n);
            end
            if(isempty(Em))
                Em=eye(n_y);
            end
            if(isempty(input_norm))
                input_norm=zeros(n_i,1)+inf;
            end
            if(isempty(state_norm))
                state_norm=zeros(n,1)+inf;
            end
            
            % Covert a scalar to a vector having the same entries.
            if(length(pn_norm)==1&&n>1)
                pn_norm=ones(n,1)*pn_norm;
                warning(['Norm type of process noise is a scalar,'...
                        ' converted to a vector with identical entries.']);
            end
            if(length(mn_norm)==1&&n_y>1)
                mn_norm=ones(n_y,1)*mn_norm;
                warning(['Norm type of measurement noise is a '...
                   'scalar, converted to a vector with identical entries.']);
            end
            if(length(f)==1&&(n+n_mode>2))
                f=ones(n,n_mode)*f;
                warning(['Additive constant for outputs is a '...
                   'scalar, converted to a matrix with identical entries.']);
            end
            if(length(g)==1&&(n_y+n_mode>2))
                g=ones(n_y,n_mode)*g;
                warning(['Additive constant for states is a scalar, '...
                    'converted to a matrix with identical entries.']);
            end
            if(length(input_norm)==1&&n_i>1)
                input_norm=ones(n_i,1)*input_norm;
                warning(['Norm type of input is a scalar, converted to'...
                    ' a vector with identical entries.']);
            end
            if(length(state_norm)==1&&n>1)
                state_norm=ones(n,1)*state_norm;
                warning(['Norm type of states is a scalar, converted'...
                    ' to a vector with identical entries.']);
            end
            
            % Check the noise parameters.
            if(length(pn_norm)~=n||~isvector(pn_norm))
                error(['The number of norm types for process noise is not'...
                      ' correct.']);
            end
            if(length(mn_norm)~=n_y||~isvector(mn_norm))
                error(['The number of norm types for measurement noise is'...
                      ' not correct.']);
            end
            
            % Check the state and input parameters.
            if(length(input_norm)~=n_i||~isvector(input_norm))
                error('The number of norm types for input is incorrect.');
            end
            if(length(state_norm)~=n||~isvector(state_norm))
                error('The number of norm types for states is incorrect.');
            end
            
            
            % Check the constant f and g.
            if(size(f,1)~=n||size(f,2)~=n_mode)
                error('The additive constant (notation f) is not valid.');
            end
            if(size(g,1)~=n_y||size(g,2)~=n_mode)
                error('The additive constant (notation g) is not valid.');
            end
            
            % Assign values after checking.
            for i=1:n_mode
                sys.mode(i).A=A(:,:,i);
                sys.mode(i).B=B(:,:,i);
                sys.mode(i).C=C(:,:,i);
                sys.mode(i).D=D(:,:,i);
                sys.mode(i).P=P(:,:,i);
                sys.mode(i).M=M(:,i);
                sys.mode(i).g=g(:,i);
                sys.mode(i).f=f(:,i);
            end
            sys.Ep=Ep;
            sys.Em=Em;
            sys.pn_norm=pn_norm;
            sys.mn_norm=mn_norm;
            sys.input_norm=input_norm;
            sys.state_norm=state_norm;
            
        end
        
    end
end
