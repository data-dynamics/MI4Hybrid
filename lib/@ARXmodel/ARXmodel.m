classdef ARXmodel
    
    % The class represents a discrete-time (possibly switched) ARX model.
    % A general discrete-time ARX model with a switching sequence sigma has
    % the following form:
    %   y[k] = (sum_i) A[sigma[k]](i)*y1[k-i] + (sum_i) C[sigma[k]](i)*u[k-i]
    %          + f[sigma[k]] + Ep*pn[k]
    %   y_n[k] = y[k] + Em*mn[k]
    % where pn is the process noise and mn is the measurement noise. A[sigma[k]]
    % and C[sigma[k]] are 3D arrays having nA and nC matrices respectively.
    %
    % Syntax:
    %   sys=ARXmodel(A,C);
    %   sys=ARXmodel(A,C,f);
    %   sys=ARXmodel(A,C,f,pn_norm,mn_norm);
    %   sys=ARXmodel(A,C,f,pn_norm,mn_norm,Ep,Em);
    %   sys=ARXmodel(A,C,f,pn_norm,mn_norm,Ep,Em,input_norm);
    %   
    % Author: MI4Hybrid
    % Date: May 22nd, 2015
    
    % Notations:
    %   n_y -- number of outputs
    %   n_mode -- number of modes
    %   n_i -- number of inputs
    
    properties(SetAccess = protected)
        % A set of discrete-time ARX modes.
        % e.g. mode(i).A and mode(i).C represent the i-th mode.
        mode
        % An n_y-by-1 column vector representing the norm types of process
        % noise.
        pn_norm
        % An n_y-by-1 column vector representing the norm types of measurement
        % noise.
        mn_norm
        % f is an n_y-by-n_mode matrix.
        f
        % Ep is an n_y-by-n_y matrix for process noise.
        Ep
        % Em is an n_y-by-n_y matrix for measurement noise.
        Em
        % An n_i-by-n_i column vector representing the norm types of inputs.
        input_norm
        % A mark (a string) stating that the model is a regular or switched
        % ARX model.
        mark
    end
    
    methods
        
        % If there is only one mode, the model is not switchable.
        function sys=ARXmodel(A,C,f,pn_norm,mn_norm,Ep,Em,input_norm)
            
            % Check A and C.
            if(size(A,4)~=size(C,4))
                error(['The first two arguments must represent the '...
                    'same number of modes.']);
            elseif(size(A,1)~=size(C,1))
                error(['The matrices described by the first two '...
                    'arguments are not consistent.']);
            elseif(size(A,1)~=size(A,2))
                error('The first argument must be square matrices.');
            else
                n_mode=size(A,4); % number of modes
                n_y=size(A,1); % number of outputs
                n_i=size(C,2); % number of inputs
            end
            
            % Indicate the model type.
            if(n_mode>1)
                sys.mark='swarx';
            elseif(n_mode==1)
                sys.mark='arx';
            else
                error('There should be at least one mode.');
            end
            
            % Set up default values if the parameters are not specified.
            if(nargin==2)
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                f=zeros(n_y,n_mode);
                Ep=eye(n_y);
                Em=eye(n_y);
                input_norm=zeros(n_i,1)+inf;
            end
            if(nargin==3)
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n_y);
                Em=eye(n_y);
                input_norm=zeros(n_i,1)+inf;
            end
            if(nargin==5)
                Ep=eye(n_y);
                Em=eye(n_y);
                input_norm=zeros(n_i,1)+inf;
            end
            if(nargin==7)
                input_norm=zeros(n_i,1)+inf;
            end
            
            % Set up default values for empty inputs.
            if(isempty(pn_norm))
                pn_norm=zeros(n_y,1)+inf;
            end
            if(isempty(mn_norm))
                mn_norm=zeros(n_y,1)+inf;
            end
            if(isempty(Ep))
                Ep=eye(n_y);
            end
            if(isempty(Em))
                Em=eye(n_y);
            end
            if(isempty(input_norm))
                input_norm=zeros(n_i,1)+inf;
            end
            
            % Covert a scalar (pn_norm or mn_norm) to a vector having the
            % same entries.
            if(length(pn_norm)==1&&n_y>1)
                pn_norm=ones(size(A,1),1)*pn_norm;
                warning(['Input norm type of process noise is a scalar,'...
                        ' converted to a vector with identical entries.']);
            end
            if(length(mn_norm)==1&&n_y>1)
                mn_norm=ones(size(A,1),1)*mn_norm;
                warning(['Input norm type of measurement noise is a '...
                    'scalar, converted to a vector with identical entries.']);
            end
            if(length(f)==1&&(n_y+n_mode>2))
                f=ones(size(A,1),n_mode)*f;
                warning(['Input additive constant for outputs is a scalar'...
                    ', converted to a matrix with identical entries.']);
            end
            if(length(input_norm)==1&&n_i>1)
                input_norm=ones(n_i,1)*input_norm;
                warning(['Input norm type is a scalar, converted to a'...
                    ' vector with identical entries.']);
            end
            
            % Check the noise parameters.
            if(length(pn_norm)~=n_y||~isvector(pn_norm))
                error(['The number of norm types for process noise is'...
                    ' not correct.']);
            elseif(length(mn_norm)~=n_y||~isvector(mn_norm))
                error(['The number of norm types for measurement noise'...
                    ' is not correct.']);
            end
            
            % Check the input parameters.
            if(length(input_norm)~=n_i||~isvector(input_norm))
                error('The number of norm types for input is incorrect.');
            end
            
            % Check Ep and Em
            if(size(Ep,1)~=n_y||size(Ep,2)~=n_y)
                error(['The factor (matrix) for process noise is not '...
                      'valid.']);
            end
            if(size(Em,1)~=n_y||size(Em,2)~=n_y)
                error(['The factor (matrix) for measurement noise is not'...
                      ' valid.']);
            end
            
            % Make the 3rd dimension of A and C to be the same by adding
            % zeros in order to easily run the simulation.
            if(size(A,3)>size(C,3))
                current_deg=size(C,3);
                n_in=size(C,2);
                add_deg=size(A,3)-size(C,3);
                for i=1:n_mode
                    for j=1:add_deg
                        C(:,:,current_deg+j,i)=zeros(n_y,n_in);
                    end
                end
            elseif(size(A,3)<size(C,3))
                current_deg=size(A,3);
                add_deg=size(C,3)-size(A,3);
                for i=1:n_mode
                    for j=1:add_deg
                        A(:,:,current_deg+j,i)=zeros(n_y,n_y);
                    end
                end
            end
            
            % Assign values after checking.
            for i=1:n_mode
                sys.mode(i).A=A(:,:,:,i);
                sys.mode(i).C=C(:,:,:,i);
            end
            sys.Ep=Ep;
            sys.Em=Em;
            sys.pn_norm=pn_norm;
            sys.mn_norm=mn_norm;
            sys.f=f;
            sys.input_norm=input_norm;
            
        end
        
    end
    
end