classdef UnPolyModel < PolyModel
    
    % This class describes a (non-switched) polynomial model with parameter
    % uncertainty. The uncertainty is specified to its coefficient matrix,
    % that is, this matrix has the form coeffmat_true = coeffmat +
    % delta_coeffmat.
    %
    % Syntax:
    %   sys=UnPolyModel(sys,d_coeffmat);
    %   sys=UnPolyModel(degmat,coeffmat,d_coeffmat);
    %   sys=UnPolyModel(degmat,coeffmat,d_coeffmat,pn_norm,mn_norm);
    %   sys=UnPolyModel(degmat,coeffmat,d_coeffmat,pn_norm,mn_norm,Ep,Em);
    %   sys=UnPolyModel(degmat,coeffmat,d_coeffmat,pn_norm,mn_norm,Ep,...
    %                   Em,input_norm);
    %   sys=UnPolyModel(degmat,coeffmat,d_coeffmat,pn_norm,mn_norm,Ep,...
    %                   Em,input_norm,state_norm);
    %
    % Author: MI4Hybrid
    % Date: June 5th, 2015
    
    % Notations:
    %   n -- number of states/outputs/polynomials
    %   n_i -- number of inputs
    %   n_m -- number of monomials in the model
    
    properties(SetAccess = protected)
        % This matrix has the same dimension as the coefficient matrix. It
        % defines the uncertainty limits for the coefficient matrix.
        d_coeffmat
    end
    
    methods
        
        function sys=UnPolyModel(varargin)
            
            % To be able to convert a system model from superclass to
            % subclass (for model transfer).
            if(nargin==2&&isa(varargin{1},'polymodel'))
                sys0=varargin{1};
                degmat=sys0.degmat;
                coeffmat=sys0.coeffmat;
                pn_norm=sys0.pn_norm;
                mn_norm=sys0.mn_norm;
                Ep=sys0.Ep;
                Em=sys0.Em;
                input_norm=sys0.input_norm;
                state_norm=sys0.state_norm;
                delta_coeffmat=varargin{2};
                [n,n_m]=size(coeffmat);
                n_i=size(degmat,2)-n;
            end
            
            % Obtain basic model information (not for model transfer).
            if(nargin>=3)
                degmat=varargin{1};
                coeffmat=varargin{2};
                delta_coeffmat=varargin{3};
                [n,n_m]=size(coeffmat);
                n_i=size(degmat,2)-n;
            end
            
            % Set up default values if parameters are not specified (not
            % for model transfer).
            if(nargin==3)
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n,1)+inf;
                Ep=eye(n);
                Em=eye(n);
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==5)
                pn_norm=varargin{4};
                mn_norm=varargin{5};
                Ep=eye(n);
                Em=eye(n);
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==7)
                pn_norm=varargin{4};
                mn_norm=varargin{5};
                Ep=varargin{6};
                Em=varargin{7};
                input_norm=zeros(n_i,1)+inf;
                state_norm=zeros(n,1)+inf;
            elseif(nargin==8)
                pn_norm=varargin{4};
                mn_norm=varargin{5};
                Ep=varargin{6};
                Em=varargin{7};
                input_norm=varargin{8};
                state_norm=zeros(n,1)+inf;
            end
            
            % Set up default values for empty arguments.
            if(isempty(delta_coeffmat))
                delta_coeffmat=zeros(n,n_m);
            end
            
            % Convert scalars to vectors.
            if(length(delta_coeffmat)==1&&length(coeffmat)~=1)
                delta_coeffmat=zeros(n,n_m)+delta_coeffmat;
                warning(['The uncertainty constraint for the '...
                    'coefficient matrix is a scalar, converted to an '...
                    'array with the same entries.']);
            end
            
            % Check if the uncertainty is consistent with coeffmat.
            [m2,n2]=size(delta_coeffmat);
            if(n~=m2||n_m~=n2)
                error(['The uncertainty constraint for the coefficient'...
                    ' matrix should have the same size as the '...
                    'coefficient matrix.']);
            end
            
            % Call the constructor of the super class.
            sys@PolyModel(degmat,coeffmat,pn_norm,mn_norm,Ep,Em,...
                input_norm,state_norm);
            
            % Assign the uncertainty constraints.
            sys.d_coeffmat=delta_coeffmat;
            
        end
        
    end
    
end