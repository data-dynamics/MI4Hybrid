classdef UnStateSpace < StateSpace
    
    % This class represents a discrete-time state-space model with parameter
    %   uncertainty, that is, A, B, C, and D matrices of a system mode are in
    %   the form of A = A_true + delta_A, B = B_true + delta_B, C = C_true +
    %   delta_C, and D = D_true + delta_D.
    %
    % Syntax:
    %   sys=UnStateSpace(sys,d_A,d_B,d_C,d_D);
    %   sys=UnStateSpace(A,B,C,D);
    %   sys=UnStateSpace(A,B,C,D,g,f);
    %   sys=UnStateSpace(A,B,C,D,g,f,d_A,d_B,d_C,d_D);
    %   sys=UnStateSpace(A,B,C,D,g,f,d_A,d_B,d_C,d_D,pn_norm,mn_norm);
    %   sys=UnStateSpace(A,B,C,D,g,f,d_A,d_B,d_C,d_D,pn_norm,mn_norm,Ep,Em);
    %
    % Author: MI4Hybrid
    % Date: May 26th, 2015
    
    properties(SetAccess=protected)
        % A set of uncertainty limits of discrete-time state-space modes.
        % e.g. d_mode(i).A, d_mode(i).B, d_mode(i).C, d_mode(i).D represent
        %   the uncertainty limits for the i-th mode.
        d_mode
    end
    
    methods
        
        % Constructor for a (switched) state space model with uncertainty.
        function sys=UnStateSpace(varargin)
            
            % To be able to convert a system model from superclass to subclass.
            if(nargin==5&&isa(varargin{1},'StateSpace'))
                sys0=varargin{1};
                n_mode=size(sys0.mode,2);
                for i=1:n_mode
                    A(:,:,i)=sys0.mode(i).A;
                    B(:,:,i)=sys0.mode(i).B;
                    C(:,:,i)=sys0.mode(i).C;
                    D(:,:,i)=sys0.mode(i).D;
                end
                g=sys0.g;
                f=sys0.f;
                pn_norm=sys0.pn_norm;
                mn_norm=sys0.mn_norm;
                Ep=sys0.Ep;
                Em=sys0.Em;
                d_A=varargin{2};
                d_B=varargin{3};
                d_C=varargin{4};
                d_D=varargin{5};
            end
            
            % Obtain basic model information (not for model transfer).
            if(nargin>=4&&nargin~=5)
                A=varargin{1};
                B=varargin{2};
                C=varargin{3};
                D=varargin{4};
                n=size(A,1); % number of states
                n_y=size(C,1); % number of outputs
                n_mode=size(A,3); % number of modes
            end
            
            % Set up default values if parameters are not specified (not for model transfer).
            if(nargin==4)
                g=zeros(n,n_mode);
                f=zeros(n_y,n_mode);
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
                d_A=zeros(size(A));
                d_B=zeros(size(B));
                d_C=zeros(size(C));
                d_D=zeros(size(D));
            end
            if(nargin==6)
                g=varargin{5};
                f=varargin{6};
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
                d_A=zeros(size(A));
                d_B=zeros(size(B));
                d_C=zeros(size(C));
                d_D=zeros(size(D));
            end
            if(nargin==10)
                g=varargin{5};
                f=varargin{6};
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=1;
                Em=1;
                d_A=varargin{7};
                d_B=varargin{8};
                d_C=varargin{9};
                d_D=varargin{10};
            end
            if(nargin==12)
                g=varargin{5};
                f=varargin{6};
                pn_norm=varargin{11};
                mn_norm=varargin{12};
                Ep=1;
                Em=1;
                d_A=varargin{7};
                d_B=varargin{8};
                d_C=varargin{9};
                d_D=varargin{10};
            end
            if(nargin==14)
                g=varargin{5};
                f=varargin{6};
                pn_norm=varargin{11};
                mn_norm=varargin{12};
                Ep=varargin{13};
                Em=varargin{14};
                d_A=varargin{7};
                d_B=varargin{8};
                d_C=varargin{9};
                d_D=varargin{10};
            end
            
            % Convert scalars to vectors.
            if(length(d_A)==1&&length(A)~=1)
                d_A=zeros(size(A))+d_A;
                warning('The uncertainty constraints for the A matrices is a scalar, converted to an array with the same entries.');
            end
            if(length(d_B)==1&&length(B)~=1)
                d_B=zeros(size(B))+d_B;
                warning('The uncertainty constraints for the B matrices is a scalar, converted to an array with the same entries.');
            end
            if(length(d_C)==1&&length(C)~=1)
                d_C=zeros(size(C))+d_C;
                warning('The uncertainty constraints for the C matrices is a scalar, converted to an array with the same entries.');
            end
            if(length(d_D)==1&&length(D)~=1)
                d_D=zeros(size(D))+d_D;
                warning('The uncertainty constraints for the D matrices is a scalar, converted to an array with the same entries.');
            end
            
            % Check if the uncertainty is consistent with parameters.
            if(~isequal(size(A),size(d_A))||~isequal(size(B),size(d_B))||~isequal(size(C),size(d_C))||~isequal(size(D),size(d_D)))
                error('The uncertainty constraints cannot match with system parameters.');
            end
            
            % Call the constructor of the super class.
            sys@StateSpace(A,B,C,D,g,f,pn_norm,mn_norm,Ep,Em);
            
            % Assign the uncertainty constraints.
            for i=1:n_mode
                sys.d_mode(i).A=d_A(:,:,i);
                sys.d_mode(i).B=d_B(:,:,i);
                sys.d_mode(i).C=d_C(:,:,i);
                sys.d_mode(i).D=d_D(:,:,i);
            end
        
        end
        
    end
    
end