classdef UnARXmodel < ARXmodel
    
    % This class represents a (possibly switched) discrete-time ARX model
    %   with parameter uncertainty, that is, A and C matrices of a system
    %   mode are in the form of A = A_true + delta_A, C = C_true + delta_C.
    %
    % Syntax:
    %   sys=UnARXmodel(sys,d_A,d_C);
    %   sys=UnARXmodel(A,C);
    %   sys=UnARXmodel(A,C,f);
    %   sys=UnARXmodel(A,C,f,d_A,d_C);
    %   sys=UnARXmodel(A,C,f,d_A,d_C,pn_norm,mn_norm);
    %   sys=UnARXmodel(A,C,f,d_A,d_C,pn_norm,mn_norm,Ep,Em);
    %
    % Author: MI4Hybrid
    % Date: May 26th, 2015
    
    properties(SetAccess=protected)
        % A set of uncertainty limits of discrete-time ARX modes.
        % e.g. d_mode(i).A, d_mode(i).C represent the uncertainty limits for
        %   the i-th mode.
        d_mode
    end
    
    methods
        
        % Constructor for a (switched) ARX model with uncertainty.
        function sys=UnARXmodel(varargin)
            
            % To be able to convert a system model from superclass to subclass.
            if(nargin==3&&isa(varargin{1},'ARXmodel'))
                sys0=varargin{1};
                n_mode=size(sys0.mode,2);
                for i=1:n_mode
                    A(:,:,i)=sys0.mode(i).A;
                    C(:,:,i)=sys0.mode(i).C;
                end
                f=sys0.f;
                pn_norm=sys0.pn_norm;
                mn_norm=sys0.mn_norm;
                Ep=sys0.Ep;
                Em=sys0.Em;
                d_A=varargin{2};
                d_C=varargin{3};
            end
            
            % Obtain basic model information (not for model transfer).
            if(nargin>=2&&~isa(varargin{1},'ARXmodel'))
                A=varargin{1};
                C=varargin{2};
                n_y=size(C,1); % number of outputs
                n_mode=size(A,3); % number of modes
            end
            
            % Set up default values if parameters are not specified (not for model transfer).
            if(nargin==2)
                f=zeros(n_y,n_mode);
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
                d_A=zeros(size(A));
                d_C=zeros(size(C));
            elseif(nargin==3&&~isa(varargin{1},'ARXmodel'))
                f=varargin{3};
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
                d_A=zeros(size(A));
                d_C=zeros(size(C));
            elseif(nargin==5)
                f=varargin{3};
                pn_norm=zeros(n_y,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
                d_A=varargin{4};
                d_C=varargin{5};
            elseif(nargin==7)
                f=varargin{3};
                pn_norm=varargin{6};
                mn_norm=varargin{7};
                Ep=1;
                Em=1;
                d_A=varargin{4};
                d_C=varargin{5};
            elseif(nargin==9)
                f=varargin{3};
                pn_norm=varargin{6};
                mn_norm=varargin{7};
                Ep=varargin{8};
                Em=varargin{9};
                d_A=varargin{4};
                d_C=varargin{5};
            end
            
            % Convert scalars to vectors.
            if(length(d_A)==1&&length(A)~=1)
                d_A=zeros(size(A))+d_A;
                warning('The uncertainty constraints for the A matrices is a scalar, converted to an array with the same entries.');
            end
            if(length(d_C)==1&&length(C)~=1)
                d_C=zeros(size(C))+d_C;
                warning('The uncertainty constraints for the C matrices is a scalar, converted to an array with the same entries.');
            end
            
            % Check if the uncertainty is consistent with parameters.
            if(~isequal(size(A),size(d_A))||~isequal(size(C),size(d_C)))
                error('The uncertainty constraints cannot match with system parameters.');
            end
            
            % Call the constructor of the super class.
            sys@ARXmodel(A,C,f,pn_norm,mn_norm,Ep,Em);
            
            % Assign the uncertainty constraints.
            for i=1:n_mode
                sys.d_mode(i).A=d_A(:,:,i);
                sys.d_mode(i).C=d_C(:,:,i);
            end
            
        end
        
    end
    
end  