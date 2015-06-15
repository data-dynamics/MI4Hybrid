classdef ARXmodel
    
    % The class represents a discrete-time (possibly switched) ARX model.
    % A general discrete-time ARX model with a switching sequence sigma has
    % the following form:
    %   y[k] = (sum) A[sigma[k]]*y1[k-nA:k-1] + (sum) C[sigma[k]]*u[k-nC:k-1] + f[sigma[k]] + Ep*pn[k]
    %   y_n[k] = y[k] + Em*mn[k]
    % where pn is the process noise and mn is the measurement noise. A[sigma[k]]
    % and C[sigma[k]] are regressors (i.e. 3D arrays) having nA and nC degrees
    % respectively.
    %
    % Syntax:
    %   sys=ARXmodel(A,C);
    %   sys=ARXmodel(A,C,f);
    %   sys=ARXmodel(A,C,f,pn_norm,mn_norm);
    %   sys=ARXmodel(A,C,f,pn_norm,mn_norm,Ep,Em);
    %   
    % Author: MI4Hybrid
    % Date: May 22nd, 2015
    
    % Notations:
    %   n_y -- number of outputs
    %   n_mode -- number of modes
    properties(SetAccess=protected)
        % A set of discrete-time ARX modes.
        % e.g. mode(i).A and mode(i).C represent the i-th mode.
        mode
        % An n_y-by-1 column vector representing the norm types of process noise.
        pn_norm
        % An n_y-by-1 column vector representing the norm types of measurement noise.
        mn_norm
        % f is an n_y-by-n_mode matrix.
        f
        % Ep is an n_y-by-n_y matrix for process noise.
        Ep
        % Em is an n_y-by-n_y matrix for measurement noise.
        Em
        % A mark (a string) stating that the model is a regular or switched ARX model.
        mark
    end
    
    methods
        
        % If there is only one mode, the model is not switchable.
        function sys=ARXmodel(A,C,f,pn_norm,mn_norm,Ep,Em)
            
            % Check A and C.
            if(size(A,4)~=size(C,4))
                error('The first two arguments must represent the same number of modes.');
            elseif(size(A,1)~=size(C,1))
                error('The matrices described by the first two arguments are not consistent.');
            elseif(size(A,1)~=size(A,2))
                error('The first argument must be square matrices.');
            else
                n_mode=size(A,4); % number of modes
                n_y=size(A,1); % number of outputs
            end
            
            % Indicate the model type.
            if(n_mode>1)
                sys.mark='swarx';
            elseif(n_mode==1)
                sys.mark='arx';
            else
                error('There must be at least one mode.');
            end
            
            % Set up default values if the parameters are not specified.
            if(nargin==2)
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                f=zeros(n_y,n_mode);
                Ep=1;
                Em=1;
            end
            if(nargin==3)
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
            end
            if(nargin==5)
                Ep=1;
                Em=1;
            end
            
            % Covert a scalar (pn_norm or mn_norm) to a vector having the same entries.
            if(length(pn_norm)==1&&n_y>1)
                pn_norm=ones(size(A,1),1)*pn_norm;
                warning('Input norm type of process noise is a scalar, converted to a vector with identical entries.');
            end
            if(length(mn_norm)==1&&n_y>1)
                mn_norm=ones(size(A,1),1)*mn_norm;
                warning('Input norm type of measurement noise is a scalar, converted to a vector with identical entries.');
            end
            if(length(f)==1&&(n_y+n_mode>2))
                f=ones(size(A,1),n_mode)*f;
                warning('Input additive constant for outputs is a scalar, converted to a matrix with identical entries.');
            end
            
            % Check the noise parameters.
            if(length(pn_norm)~=n_y)
                error('The number of norm types for process noise is not correct.');
            elseif(length(mn_norm)~=n_y)
                error('The number of norm types for measurement noise is not correct.');
            end
            
            % Check Ep and Em
            if(Ep~=1&&(size(Ep,1)~=n_y||size(Ep,2)~=n_y))
                error('The factor (matrix) for process noise is not correct.');
            end
            if(Em~=1&&(size(Em,1)~=n_y||size(Em,2)~=n_y))
                error('The factor (matrix) for measurement noise is not correct.');
            end
            
            % Make the 3rd dimension of A and C to be the same by adding zero
            % columns in order to easily run the simulation.
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
            
        end
        
    end
    
end