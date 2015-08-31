classdef PolyModel
    
    % This class describes a discrete-time (not switched) polynomial model.
    % A general polynomial model has the form:
    %   x[k+1] = P(x[k],u[k]) + Ep*pn[k]
    %   y[k] = x[k] + Em*mn[k]
    % where pn is the process noise and mn is the measurement noise.
    %
    % Syntax:
    %   sys=PolyModel(degmat,coeffmat);
    %   sys=PolyModel(degmat,coeffmat,pn_norm,mn_norm);
    %   sys=PolyModel(degmat,coeffmat,pn_norm,mn_norm,Ep,Em);
    %   sys=PolyModel(degmat,coeffmat,pn_norm,mn_norm,Ep,Em,input_norm);
    %   sys=PolyModel(degmat,coeffmat,pn_norm,mn_norm,Ep,Em,input_norm,...
    %                   state_norm);
    % Each row of degmat represents a monomial. The number of coulumns of
    % degmat is the number of states plus the number of inputs. Each row of
    % coeffmat represent a polynomial. The number of columns of coeffmat is
    % the number of monomials.
    %
    % Author: MI4Hybrid
    % Date: May 29th, 2015
    
    % Notations:
    %   n -- number of states/outputs/polynomials
    %   n_i -- number of inputs
    %   n_m -- number of monomials in the model
    
    properties(SetAccess = protected)
        % An n_m-by-(n_i+n) matrix.
        degmat
        % An n-by-n_m matrix.
        coeffmat
        % An n-by-1 column vector representing the norm types of process
        % noise.
        pn_norm
        % An n-by-1 column vector representing the norm types of 
        % measurement noise.
        mn_norm
        % Ep is an n-by-n matrix.
        Ep
        % Em is an n-by-n matrix.
        Em
        % An n_i-by-1 column vector representing the norm types of input.
        input_norm
        % An n-by-1 column vector representing the norm types of state.
        state_norm
        % A mark (a string) stating the model type.
        mark
    end
    
    methods
        
        function sys=PolyModel(degmat,coeffmat,pn_norm,mn_norm,Ep,Em,...
                input_norm,state_norm)
            
            % Check the consistency of degmat and coeffmat.
            [n,n_m]=size(coeffmat);
            n_i=size(degmat,2)-n;
            if(size(degmat,1)~=n_m||n_i<0)
                error(['The degree matrix and coefficient matrix are '...
                      'not consistent.']);
            else
                sys.mark='poly';
            end
            
            % Set up default values if parameters are not specified.
            if(nargin==2)
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n,1)+inf;
                Ep=eye(n);
                Em=eye(n);
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==4)
                Ep=eye(n);
                Em=eye(n);
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==6)
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==7)
                state_norm=zeros(n,1)+inf;
            end
            
            % Set up default values for empty inputs.
            if(isempty(pn_norm))
                pn_norm=zeros(n,1)+inf;
            end
            if(isempty(mn_norm))
                mn_norm=zeros(n,1)+inf;
            end
            if(isempty(Ep))
                Ep=eye(n);
            end
            if(isempty(Em))
                Em=eye(n);
            end
            if(isempty(input_norm))
                input_norm=zeros(n_i,1)+inf;
            end
            if(isempty(state_norm))
                state_norm=zeros(n,1)+inf;
            end
            
            % Convert a scalar (pn_norm or mn_norm) to a vector having the
            % same entries.
            if(length(pn_norm)==1&&n>1)
                pn_norm=ones(n,1)*pn_norm;
                warning(['Norm type of process noise is a scalar, '...
                    'converted to a vector with identical entries.']);
            end
            if(length(mn_norm)==1&&n>1)
                mn_norm=ones(n,1)*mn_norm;
                warning(['Norm type of measurement noise is a scalar, '...
                    'converted to a vector with identical entries.']);
            end
            if(length(input_norm)==1&&n_i>1)
                input_norm=ones(n_i,1)*input_norm;
                warning(['Norm type of input is a scalar, '...
                    'converted to a vector with identical entries.']);
            end
            if(length(state_norm)==1&&n>1)
                state_norm=ones(n,1)*state_norm;
                warning(['Norm type of state is a scalar, '...
                    'converted to a vector with identical entries.']);
            end
            
            % Check the noise parameters.
            if(length(pn_norm)~=n||~isvector(pn_norm))
                error(['The number of norm types for process noise is'...
                      ' incorrect.']);
            end
            if(length(mn_norm)~=n||~isvector(mn_norm))
                error(['The number of norm types for measurement noise'...
                      ' is incorrect.']);
            end
            
            % Check the input parameters.
            if(length(input_norm)~=n_i||~isvector(input_norm))
                error('The number of norm types for input is incorrect.');
            end
            
            % Chenk the state parameters.
            if(length(state_norm)~=n||~isvector(state_norm))
                error('The number of norm types for state is incorrect.');
            end
            
            % Check Ep and Em
            if(size(Ep,1)~=n||size(Ep,2)~=n)
                error(['The factor (matrix) for process noise is not '...
                    'valid.']);
            end
            if(size(Em,1)~=n||size(Em,2)~=n)
                error(['The factor (matrix) for measurement noise is '...
                    'not valid.']);
            end
            
            % Assign values to creat the model.
            sys.degmat=degmat;
            sys.coeffmat=coeffmat;
            sys.pn_norm=pn_norm;
            sys.mn_norm=mn_norm;
            sys.Ep=Ep;
            sys.Em=Em;
            sys.input_norm=input_norm;
            sys.state_norm=state_norm;
            
        end
        
    end
    
end