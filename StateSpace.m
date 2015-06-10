classdef StateSpace
    
    % The class represents a discrete-time (possibly switched) state-space model.
    % A general discrete-time state-space model has the form:
    %   x[k+1] = A[k]*x[k] + B[k]*u[k] + g + Ep*pn[k]
    %   y[k] = C[k]*x[k] + D[k]*u[k] + f
    %   y_n[k] = y[k] + Em*mn[k]
    % where pn is the process noise and mn is the measurement noise.
    %
    % Constructor syntax:
    %   sys=StateSpace(A,B,C,D);
    %   sys=StateSpace(A,B,C,D,g,f);
    %   sys=StateSpace(A,B,C,D,g,f,pn_norm,mn_norm);
    %   sys=StateSpace(A,B,C,D,g,f,pn_norm,mn_norm,Ep,Em);
    %
    % Author: MI4Hybrid
    % Date: May 22nd, 2015

    % Notations:
    %   n -- number of states
    %   n_y -- number of outputs
    %   n_mode -- number of modes
    properties(SetAccess=protected)
        % A set of discrete-time state-space modes.
        % e.g. mode(i).A, mode(i).B, mode(i).C, and mode(i).D represent the i-th mode.
        mode
        % An n-by-1 column vector representing the norm types of process noise.
        pn_norm
        % An n_y-by-1 column vector representing the norm types of measurement noise.
        mn_norm
        % Ep is an n-by-n matrix.
        Ep
        % Em is an n_y-by-n_y matrix.
        Em
        % g is an n-by-n_mode matrix.
        g
        % f is an n_y-by-n_mode matrix.
        f
        % A mark (a string) stating that the model is a regular or switched state-space model.
        mark
    end
    
    methods
        
        % If there is only one mode, the model is not switchable.
        function sys=StateSpace(A,B,C,D,g,f,pn_norm,mn_norm,Ep,Em)
            
            % Check A, B, C, and D.
            if(size(A,3)~=size(B,3)||size(A,3)~=size(C,3)||size(A,3)~=size(D,3))
                error('A, B, C, and D must have the same number of matrices.');
            elseif(size(A,1)~=size(A,2)||size(A,1)~=size(B,1))
                error('A must be a (or a set of) square matrices, or A and B not consistent.');
            elseif(size(C,2)~=size(A,1)||size(C,1)~=size(D,1))
                error('C is not consistent with A or D.');
            else
                n_mode=size(A,3); % number of modes
                n_y=size(C,1); % number of outputs
                n=size(A,1); % number of states
            end
            
            % Indicate the model type.
            if(n_mode>1)
                sys.mark='swss';
            elseif(n_mode==1)
                sys.mark='ss';
            else
                error('There must be at least one mode.');
            end
            
            % Set up default values if parameters are not specified.
            if(nargin==4)
                g=zeros(n,n_mode);
                f=zeros(n_y,n_mode);
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
            end
            if(nargin==6)
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
            end
            if(nargin==8)
                Ep=1;
                Em=1;
            end
            
            % Covert a scalar (pn_norm or mn_norm) to a vector having the same entries.
            if(length(pn_norm)==1&&n>1)
                pn_norm=ones(n,1)*pn_norm;
                warning('Input norm type of process noise is a scalar, converted to a vector with identical entries.');
            end
            if(length(mn_norm)==1&&n_y>1)
                mn_norm=ones(n_y,1)*mn_norm;
                warning('Input norm type of measurement noise is a scalar, converted to a vector with identical entries.');
            end
            if(length(f)==1&&(n_y+n_mode>2))
                f=ones(n_y,n_mode)*f;
                warning('Input additive constant for outputs is a scalar, converted to a matrix with identical entries.');
            end
            if(length(g)==1&&(n+n_mode>2))
                g=ones(n,n_mode)*g;
                warning('Input additive constant for states is a scalar, converted to a matrix with identical entries.');
            end
            
            % Check the noise parameters.
            if(length(pn_norm)~=n)
                error('The number of norm types for process noise is not correct.');
            end
            if(length(mn_norm)~=n_y)
                error('The number of norm types for measurement noise is not correct.');
            end
            
            % Check Ep and Em
            if(Ep~=1&&(size(Ep,1)~=n||size(Ep,2)~=n))
                error('The factor (matrix) for process noise is not correct.');
            end
            if(Em~=1&&(size(Em,1)~=n_y||size(Em,2)~=n_y))
                error('The factor (matrix) for measurement noise is not correct.');
            end
            
            % Assign values after checking.
            for i=1:n_mode
                sys.mode(i).A=A(:,:,i);
                sys.mode(i).B=B(:,:,i);
                sys.mode(i).C=C(:,:,i);
                sys.mode(i).D=D(:,:,i);
            end
            sys.Ep=Ep;
            sys.Em=Em;
            sys.pn_norm=pn_norm;
            sys.mn_norm=mn_norm;
            sys.g=g;
            sys.f=f;
            
        end
        
    end
    
end