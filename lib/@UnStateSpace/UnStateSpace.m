classdef UnStateSpace < StateSpace
    
    % This class represents a (possibly switched) discrete-time state-space
    % model with parameter uncertainty. That is, each system mode has its
    % parameters in the form of A_true = A + delta_A, B_true = B + delta_B,
    % C_true = C + delta_C, D_true = D + delta_D, g_true = g + delta_g, and
    % f_true = f + delta_f.
    %
    % Syntax:
    %   sys=UnStateSpace(sys,d_A,d_B,d_C,d_D,d_g,d_f);
    %   sys=UnStateSpace(A,B,C,D);
    %   sys=UnStateSpace(A,B,C,D,g,f);
    %   sys=UnStateSpace(A,B,C,D,g,f,d_A,d_B,d_C,d_D,d_g,d_f);
    %   sys=UnStateSpace(A,B,C,D,g,f,d_A,d_B,d_C,d_D,d_g,d_f,pn_norm,...
    %                    mn_norm);
    %   sys=UnStateSpace(A,B,C,D,g,f,d_A,d_B,d_C,d_D,d_g,d_f,pn_norm,...
    %                    mn_norm,Ep,Em);
    %   sys=UnStateSpace(A,B,C,D,g,f,d_A,d_B,d_C,d_D,d_g,d_f,pn_norm,...
    %                    mn_norm,Ep,Em,input_norm);
    %   sys=UnStateSpace(A,B,C,D,g,f,d_A,d_B,d_C,d_D,d_g,d_f,pn_norm,...
    %                    mn_norm,Ep,Em,input_norm,state_norm);
    %
    % Author: Z. Luo, F. Harirchi and N. Ozay
    % Date: May 26th, 2015
    
    properties(SetAccess = protected)
        % A set of uncertainty limits of discrete-time state-space modes.
        % e.g. d_mode(i).A, d_mode(i).B, d_mode(i).C, d_mode(i).D,
        % d_mode(i).g, and d_mode(i).f represent the uncertainty limits for
        % the i-th mode.
        d_mode
    end
    
    methods
        
        % Constructor for a (switched) state space model with uncertainty.
        function sys=UnStateSpace(varargin)
            
            % To be able to convert a system model from superclass to
            % subclass.
            if(nargin==7&&isa(varargin{1},'StateSpace'))
                sys0=varargin{1};
                n_mode=size(sys0.mode,2);
                for i=1:n_mode
                    A(:,:,i)=sys0.mode(i).A;
                    B(:,:,i)=sys0.mode(i).B;
                    C(:,:,i)=sys0.mode(i).C;
                    D(:,:,i)=sys0.mode(i).D;
                    g(:,i)=sys0.mode(i).g; % each column for each mode
                    f(:,i)=sys0.mode(i).f; % each column for each mode
                end
                pn_norm=sys0.pn_norm;
                mn_norm=sys0.mn_norm;
                Ep=sys0.Ep;
                Em=sys0.Em;
                d_A=varargin{2};
                d_B=varargin{3};
                d_C=varargin{4};
                d_D=varargin{5};
                d_g=varargin{6};
                d_f=varargin{7};
                n=size(A,1); % number of states
                n_y=size(C,1); % number of outputs
                n_i=size(B,2); % number of inputs
            end
            
            % Obtain basic model information (not for model transfer).
            if(nargin>=4&&~isa(varargin{1},'StateSpace'))
                A=varargin{1};
                B=varargin{2};
                C=varargin{3};
                D=varargin{4};
                n=size(A,1); % number of states
                n_y=size(C,1); % number of outputs
                n_mode=size(A,3); % number of modes
                n_i=size(B,2); % number of inputs
            end
            
            % Set up default values if parameters are not specified (not
            % for model transfer).
            if(nargin==4)
                g=zeros(n,n_mode);
                f=zeros(n_y,n_mode);
                d_A=zeros(n,n,n_mode);
                d_B=zeros(n,n_i,n_mode);
                d_C=zeros(n_y,n,n_mode);
                d_D=zeros(n_y,n_i,n_mode);
                d_g=zeros(n,n_mode);
                d_f=zeros(n_y,n_mode);
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n);
                Em=eye(n_y);
                state_norm=zeros(n,1)+inf;
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==6)
                g=varargin{5};
                f=varargin{6};
                d_A=zeros(n,n,n_mode);
                d_B=zeros(n,n_i,n_mode);
                d_C=zeros(n_y,n,n_mode);
                d_D=zeros(n_y,n_i,n_mode);
                d_g=zeros(n,n_mode);
                d_f=zeros(n_y,n_mode);
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n);
                Em=eye(n_y);
                state_norm=zeros(n,1)+inf;
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==12)
                g=varargin{5};
                f=varargin{6};
                d_A=varargin{7};
                d_B=varargin{8};
                d_C=varargin{9};
                d_D=varargin{10};
                d_g=varargin{11};
                d_f=varargin{12};
                pn_norm=zeros(n,1)+inf;
                mn_norm=zeros(n_y,1)+inf;
                Ep=eye(n);
                Em=eye(n_y);
                state_norm=zeros(n,1)+inf;
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==14)
                g=varargin{5};
                f=varargin{6};
                d_A=varargin{7};
                d_B=varargin{8};
                d_C=varargin{9};
                d_D=varargin{10};
                d_g=varargin{11};
                d_f=varargin{12};
                pn_norm=varargin{13};
                mn_norm=varargin{14};
                Ep=eye(n);
                Em=eye(n_y);
                state_norm=zeros(n,1)+inf;
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==16)
                g=varargin{5};
                f=varargin{6};
                d_A=varargin{7};
                d_B=varargin{8};
                d_C=varargin{9};
                d_D=varargin{10};
                d_g=varargin{11};
                d_f=varargin{12};
                pn_norm=varargin{13};
                mn_norm=varargin{14};
                Ep=varargin{15};
                Em=varargin{16};
                state_norm=zeros(n,1)+inf;
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==17)
                g=varargin{5};
                f=varargin{6};
                d_A=varargin{7};
                d_B=varargin{8};
                d_C=varargin{9};
                d_D=varargin{10};
                d_g=varargin{11};
                d_f=varargin{12};
                pn_norm=varargin{13};
                mn_norm=varargin{14};
                Ep=varargin{15};
                Em=varargin{16};
                state_norm=varargin{17};
                input_norm=zeros(n_i,1)+inf;
            elseif(nargin==18)
                g=varargin{5};
                f=varargin{6};
                d_A=varargin{7};
                d_B=varargin{8};
                d_C=varargin{9};
                d_D=varargin{10};
                d_g=varargin{11};
                d_f=varargin{12};
                pn_norm=varargin{13};
                mn_norm=varargin{14};
                Ep=varargin{15};
                Em=varargin{16};
                state_norm=varargin{17};
                input_norm=varargin{18};
            end
            
            % Set up default values for empty inputs (for either model
            % transfer or non-transfer).
            if(isempty(d_A))
                d_A=zeros(n,n,n_mode);
            end
            if(isempty(d_B))
                d_B=zeros(n,n_i,n_mode);
            end
            if(isempty(d_C))
                d_C=zeros(n_y,n,n_mode);
            end
            if(isempty(d_D))
                d_D=zeros(n_y,n_i,n_mode);
            end
            if(isempty(d_g))
                d_g=zeros(n,n_mode);
            end
            if(isempty(d_f))
                d_f=zeros(n_y,n_mode);
            end
            
            % Convert scalars to vectors.
            if(length(d_A)==1&&length(A)~=1)
                d_A=zeros(size(A))+d_A;
                warning(['The uncertainty constraints for the property'...
                    ' mode(i).A is a scalar, converted to an array with'...
                    ' the same entries.']);
            end
            if(length(d_B)==1&&length(B)~=1)
                d_B=zeros(size(B))+d_B;
                warning(['The uncertainty constraints for the property'...
                    ' mode(i).B is a scalar, converted to an array with'...
                    ' the same entries.']);
            end
            if(length(d_C)==1&&length(C)~=1)
                d_C=zeros(size(C))+d_C;
                warning(['The uncertainty constraints for the property'...
                    ' mode(i).C is a scalar, converted to an array with'...
                    ' the same entries.']);
            end
            if(length(d_D)==1&&length(D)~=1)
                d_D=zeros(size(D))+d_D;
                warning(['The uncertainty constraints for the property'...
                    ' mode(i).D is a scalar, converted to an array with'...
                    ' the same entries.']);
            end
            if(length(d_g)==1&&n>1)
                d_g=zeros(n,n_mode)+d_g;
                warning(['The uncertainty constraints for the additive'...
                       ' constant (for property mode(i).g) is a scalar,'...
                       'converted to an array with the same entries.']);
            end
            if(length(d_f)==1&&n_y>1)
                d_f=zeros(n_y,n_mode)+d_f;
                warning(['The uncertainty constraints for the additive'...
                       ' constant (for property mode(i).f) is a scalar,'...
                       'converted to an array with the same entries.']);
            end
            
            % Check if the uncertainty is consistent with parameters.
            if(~isequal(size(A),size(d_A))||~isequal(size(B),size(d_B))||...
                ~isequal(size(C),size(d_C))||~isequal(size(D),size(d_D)))
                error(['The uncertainty constraints (for property '...
                    'mode(i).A, mode(i).B, mode(i).C, or mode(i).D) '...
                    'cannot match with system parameters.']);
            end
            if(length(d_g)~=n||~isvector(d_g))
                error(['The uncertainty constraints (for property '...
                    'mode(i).g) cannot match with system parameters.']);
            end
            if(length(d_f)~=n_y||~isvector(d_f))
                error(['The uncertainty constraints (for property '...
                    'mode(i).f) cannot match with system parameters.']);
            end
            
            % Call the constructor of the superclass.
            sys@StateSpace(A,B,C,D,g,f,pn_norm,mn_norm,Ep,Em,input_norm,...
                state_norm);
            
            % Assign the uncertainty constraints.
            for i=1:n_mode
                sys.d_mode(i).A=d_A(:,:,i);
                sys.d_mode(i).B=d_B(:,:,i);
                sys.d_mode(i).C=d_C(:,:,i);
                sys.d_mode(i).D=d_D(:,:,i);
                sys.d_mode(i).g=d_g(:,i);
                sys.d_mode(i).f=d_f(:,i);
            end
        
        end
        
    end
    
end