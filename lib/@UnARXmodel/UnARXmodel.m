classdef UnARXmodel < ARXmodel
    
    % This class represents a (possibly switched) discrete-time ARX model
    % with parameter uncertainty. That is, A and C matrices of each system
    % mode are in the form of A_true = A + delta_A, C_true = C + delta_C.
    % In addition, f (a vector for each mode/submodel) is in the form of
    % f_true = f + delta_f.
    %
    % Syntax:
    %   sys=UnARXmodel(sys,d_A,d_C,d_f);
    %   sys=UnARXmodel(A,C);
    %   sys=UnARXmodel(A,C,f);
    %   sys=UnARXmodel(A,C,f,d_A,d_C,d_f);
    %   sys=UnARXmodel(A,C,f,d_A,d_C,d_f,pn_norm,mn_norm);
    %   sys=UnARXmodel(A,C,f,d_A,d_C,d_f,pn_norm,mn_norm,Ep,Em);
    %   sys=UnARXmodel(A,C,f,d_A,d_C,d_f,pn_norm,mn_norm,Ep,Em,input_norm);
    %
    % Author: Z. Luo, F. Harirchi and N. Ozay
    % Date: May 26th, 2015
    
    properties(SetAccess = protected)
        % A set of uncertainty limits of discrete-time ARX modes.
        % e.g. d_mode(i).A, d_mode(i).C, and d_mode(i).f represent the 
        % uncertainty limits for the i-th mode.
        d_mode
    end
    
    methods
        
        % Constructor for a (switched) ARX model with uncertainty.
        function sys=UnARXmodel(varargin)
            
            % To be able to convert a system model from superclass to 
            % subclass.
            if(nargin==4&&isa(varargin{1},'ARXmodel'))
                sys0=varargin{1};
                n_mode=size(sys0.mode,2);
                for i=1:n_mode
                    A(:,:,:,i)=sys0.mode(i).A;
                    C(:,:,:,i)=sys0.mode(i).C;
                    f(:,i)=sys0.mode(i).f; % each column for each mode
                end
                pn_norm=sys0.pn_norm;
                mn_norm=sys0.mn_norm;
                Ep=sys0.Ep;
                Em=sys0.Em;
                d_A=varargin{2};
                d_C=varargin{3};
                d_f=varargin{4};
                n_y=size(C,1); % number of outputs
                n_i=size(C,2); % number of inputs
            end
            
            % Obtain basic model information (not for model transfer).
            if(nargin>=2&&~isa(varargin{1},'ARXmodel'))
                A=varargin{1};
                C=varargin{2};
                n_y=size(C,1); % number of outputs
                n_mode=size(A,4); % number of modes
                n_i=size(C,2); % number of inputs
            end
            
            % Set up default values if parameters are not specified (not
            % for model transfer).
            if(nargin==2)
                f=zeros(n_y,n_mode);
                d_A=zeros(size(A));
                d_C=zeros(size(C));
                d_f=zeros(n_y,n_mode);
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n_y);
                Em=eys(n_y);
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==3)
                f=varargin{3};
                d_A=zeros(size(A));
                d_C=zeros(size(C));
                d_f=zeros(n_y,n_mode);
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n_y);
                Em=eys(n_y);
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==6)
                f=varargin{3};
                d_A=varargin{4};
                d_C=varargin{5};
                d_f=varargin{6};
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n_y);
                Em=eys(n_y);
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==8)
                f=varargin{3};
                d_A=varargin{4};
                d_C=varargin{5};
                d_f=varargin{6};
                pn_norm=varargin{7};
                mn_norm=varargin{8};
                Ep=eye(n_y);
                Em=eys(n_y);
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==10)
                f=varargin{3};
                d_A=varargin{4};
                d_C=varargin{5};
                d_f=varargin{6};
                pn_norm=varargin{7};
                mn_norm=varargin{8};
                Ep=varargin{9};
                Em=varargin{10};
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==11)
                f=varargin{3};
                d_A=varargin{4};
                d_C=varargin{5};
                d_f=varargin{6};
                pn_norm=varargin{7};
                mn_norm=varargin{8};
                Ep=varargin{9};
                Em=varargin{10};
                input_norm=varargin{11};
            end
            
            % Set up default values for empty inputs (for either model
            % transfer or non-transfer).
            if(isempty(d_A))
                d_A=zeros(size(A));
            end
            if(isempty(d_C))
                d_C=zeros(size(C));
            end
            if(isempty(d_f))
                d_f=zeros(n_y,n_mode);
            end
            
            % Convert scalars to matrices/vectors.
            if(length(d_A)==1&&length(A)~=1)
                d_A=zeros(size(A))+d_A;
                warning(['The uncertainty constraints for property '...
                    'model(i).A is a scalar, converted to an array '...
                    'with the same entries.']);
            end
            if(length(d_C)==1&&length(C)~=1)
                d_C=zeros(size(C))+d_C;
                warning(['The uncertainty constraints for property '...
                    'model(i).C is a scalar, converted to an array '...
                    'with the same entries.']);
            end
            if(length(d_f)==1&&n_y>1)
                d_f=zeros(n_y,n_mode)+d_f;
                warning(['The uncertainty constraints for the additive'...
                    ' constant (for property mode(i).f) is a scalar, '...
                    'converted to an array with the same entries.']);
            end
            
            % Check if the uncertainty is consistent with parameters.
            if(~isequal(size(A),size(d_A))||~isequal(size(C),size(d_C)))
                error(['The uncertainty constraints (for property '...
                    'mode(i).A or mode(i).C) cannot match with system'...
                    ' parameters.']);
            end
            if(length(d_f)~=n_y||~isvector(d_f))
                error(['The uncertainty constraints (for property '...
                    'mode(i).f) cannot match with system parameters.']);
            end
            
            % Call the constructor of the super class.
            sys@ARXmodel(A,C,f,pn_norm,mn_norm,Ep,Em,input_norm);
            
            % Assign the uncertainty constraints.
            for i=1:n_mode
                sys.d_mode(i).A=d_A(:,:,:,i);
                sys.d_mode(i).C=d_C(:,:,:,i);
                sys.d_mode(i).f=d_f(:,i);
            end
            
        end
        
    end
    
end