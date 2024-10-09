classdef DynamicalSystem < matlab.mixin.SetGetExactNames
    % DynamicalSystem construct a dynamical system object in first order or
    % second order form
    
    % M\ddot{x} + C\dot{x} + K x + fnl(x,xd) = fext(t) - Second order
    % B \dot{z} = F(z) + Fext(t)                    - First order
    % Here fnl(x) is a polynomial function of degree two or higher, which
    % is stored as a cell array such that fnl{k} corresponds to polynomials
    % of degree k+1. fnl{k} is given by a tensor of order k+2, where the
    % first mode corresponds to indices for the force vector.
    % Likewise, F(z) is a polynomial function of degree one or higher,
    % i.e., F(z) = Az + Higher order terms. F is stored as a cell array,
    % where the i-th entry gives the tensor/multiindex representation of
    % polynomials of degree i.
    
    % The second order form is converted into the first order form with z =
    % [x;\dot{x}], B = [C M;M 0], A = [-K 0;0 M], F(z)=Az+[-fnl(x,xd);0], and
    % Fext(t) = [fext(t); 0]
    
    properties
        M = []
        C = []
        K = []
        A = []
        B = []
        BinvA
        
        %% Internal Forces
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Intrusive Nonlinearity
        fnl = []            % second order nonlinear internal forces
        dfnl= []            % jacobian of internal forces
        F   = []            % first order nonlinear internal forces
        dF  = []            % jacobian of internal forces
                
        %   Semi-intrusive Nonlinearity        
        fnl_semi    = []    % function handles and parameters to compute nonlinearity, has to take N dim vectors as input
        dfnl_semi   = []    % Handle for nonlinearity derivative
                
        F_semi  = []        % function handles for first order dynamical system
        dF_semi = []        % function handle for nonlinearity derivative
        F_semi_sym  = false % by default an asymmetric function handle is assumed        
        
        %   Non-intrusive Nonlinearity
        fnl_non = [];       % second order internal nonlinear force function                     
        dfnl_non  = [];     % second order internal nonlinear force jacobian function
        
        F_non   = [];       % first order internal nonlinear force function
        dF_non  = [];       % first order internal nonlinear force jacobian function
        
        nl_input_dim     = [];     % size of input vector for nonlinearity
        
        %% Other properties
        fext = []
        Fext = []
        Omega = []
        
        n                   % dimension for x
        N                   % dimension for z
        order = [];         % whether second-order or first-order system - mandatory property
        degree              % degree of (polynomial) nonlinearity of the rhs of the dynamical system
        nKappa              % Fourier Series expansion order for Fext
        kappas =   []       % matrix with all kappas in its rows
        
        spectrum = []       % data structure constructed by linear_spectral_analysis method
        Options = DSOptions()

    end
    
    methods

        %% Constructor function

        function obj = DynamicalSystem(order)
            if nargin < 1
                error('Please input the order of the Dynamical System');
            end
            obj.order = order;
        end

        %% SET methods - general

        function set.A(obj,A)
            obj.A = A;
        end        
                  
        %% GET methods - general
              
        function A     = get.A(obj)
            A     = get_A(obj,obj.A);
        end
        
        function B     = get.B(obj)
            B     = get_B(obj,obj.B);
        end
        
        function BinvA = get.BinvA(obj)
            BinvA = get_BinvA(obj,obj.BinvA);
        end
                    
        function n = get.n(obj)
            n = length(obj.M);
        end
        
        function N     = get.N(obj)
            N = length(obj.A);
        end

        
            
        %%  Intrusive nonlinearity

        % SET Methods
        function set.fnl(obj,fnl)
            obj.fnl = set_fnl(obj,fnl);
        end

        function set.dfnl(obj,dfnl)
            obj.dfnl = set_dfnl(obj,dfnl);
        end

        function set.F(obj,F)
            obj.F = set_F(obj,F);
        end

        function set.dF(obj,dF)
            obj. dF = set_dF(obj,dF);
        end

        % GET Methods
        function F = get.F(obj)
            F = get_F(obj,obj.F);
        end
        
        function dF = get.dF(obj)
            dF = get_dF(obj,obj.dF);
        end

        function dfnl = get.dfnl(obj)
            dfnl = get_dfnl(obj,obj.dfnl);
        end

        %% Semi-intrusive nonlinearity

        % SET Methods
        function set.fnl_semi(obj,fnl_semi)
                obj.fnl_semi = fnl_semi;
        end

        function set.dfnl_semi(obj,dfnl_semi)
            obj.dfnl_semi = dfnl_semi;
        end

        function set.F_semi(obj,F_semi)
            obj.F_semi = set_F_semi(obj,F_semi);
        end
        
        function set.dF_semi(obj,dF_semi)
            obj.dF_semi = set_dF_semi(obj,dF_semi);
        end

        % GET Methods
        function fnl_semi  = get.fnl_semi(obj)
            fnl_semi = obj.fnl_semi;
        end
        
        function dfnl_semi = get.dfnl_semi(obj)
            dfnl_semi = obj.dfnl_semi;
        end
        
        function F_semi    = get.F_semi(obj)
            F_semi    = get_F_semi(obj,obj.F_semi);
        end
        
        function dF_semi   = get.dF_semi(obj)
            dF_semi   = get_dF_semi(obj,obj.dF_semi);
        end
        
        %%  Non-intrusive nonlinearity

        % SET Methods
        function set.fnl_non(obj,fnl_non)
            obj.fnl_non = fnl_non;
        end
        
        function set.dfnl_non(obj,dfnl_non)
            obj.dfnl_non = dfnl_non;
        end

        function set.F_non(obj,F_non)
            obj.F_non = F_non;
        end
        
        function set.dF_non(obj,dF_non)
            obj.dF_non = dF_non;
        end


        % GET Methods
        function fnl_non  = get.fnl_non(obj)
            fnl_non = obj.fnl_non;
        end        

        function dfnl_non = get.dfnl_non(obj)
            dfnl_non = obj.dfnl_non;
        end     
        
        function F_non    = get.F_non(obj)
            F_non = get_F_non(obj,obj.F_non);
        end

        function dF_non   = get.dF_non(obj)
            dF_non = get_dF_non(obj,obj.dF_non);
        end
        
        %% External forcing
        
        function nKappa = get.nKappa(obj)
            nKappa = numel(obj.Fext.data);
        end
        
        function kappas = get.kappas(obj)
            kappas = get_kappas(obj);
        end
        
        function Fext   = get.Fext(obj)
            Fext = get_Fext(obj,obj.Fext);                           
        end
               
        function degree = get.degree(obj)
            degree = get_degree(obj);
        end
        
        %% other methods
        
        function nl_input_dim = get.nl_input_dim(obj)
            nl_input_dim = get_nl_input_dim(obj,obj.nl_input_dim);
        end



        add_forcing(obj,f,varargin)
        [V, D, W]                   = linear_spectral_analysis(obj)
        fext                        = compute_fext(obj,t,x,xd)
        Fext                        = evaluate_Fext(obj,t,z)
        fnl                         = compute_fnl(obj,x,xd)
        dfnl                        = compute_dfnldx(obj,x,xd)
        dfnl                        = compute_dfnldxd(obj,x,xd)
        Fnl                         = evaluate_Fnl(obj,z)
        f                           = odefun(obj,t,z)
        [r, drdqdd,drdqd,drdq, c0]  = residual(obj, q, qd, qdd, t)

        % For checking input dimension
    end
end
