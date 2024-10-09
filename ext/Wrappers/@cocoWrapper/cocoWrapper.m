classdef cocoWrapper < matlab.mixin.SetGet
    %COCOWRAPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        system     % dynamical system object 
        nCycles    % number of periods in (forward) steady state simulation
        outdof     % out dof for result visualization
        periodsRatio = 1  % response period / excitation period
        Options = cocoOptions();
        initialGuess = 'forward' % 'linear', 'given'
        forwardSolver = 'Newmark' % 'Galpha'
        atlasAlg = '1d' % 'kd'
        branchSwitch = false % true
        multiFnl = []
        optnorm = 'l2' % linf
    end
    
    methods
        % constructor
        function obj = cocoWrapper(sys, ncyles, outdof)
            %COCOWRAPPER Construct an instance of this calss
            obj.system  = sys;
            obj.nCycles = ncyles;
            obj.outdof  = outdof;
        end
        
        function set.initialGuess(obj,initialGuess)
            switch lower(initialGuess)
                case 'forward'
                    obj.initialGuess = 'forward';
                case 'linear'
                    obj.initialGuess = 'linear';
                case 'given'
                    obj.initialGuess = 'given';
                otherwise
                    error('Unknown initial solver type: set forward or linear or given as intial solver types')
            end
        end 
        
        function set.forwardSolver(obj,forwardSolver)
            switch lower(forwardSolver)
                case 'newmark'
                    obj.forwardSolver = 'newmark';
                case 'galpha'
                    obj.forwardSolver = 'galpha';
                otherwise
                    error('Unknown forward solver type: set Newmark or Galpha as forward solver types')
            end
        end 
        
        function set.atlasAlg(obj,atlasAlg)
            switch lower(atlasAlg)
                case '1d'
                    obj.atlasAlg = '1d';
                case 'kd'
                    obj.atlasAlg = 'kd';
                otherwise
                    error('Unknown atlas algorithm: set 1d or kd as atlas algorithm');
            end
        end
        
        function set.optnorm(obj,optnorm)
            switch lower(optnorm)
                case 'l2'
                    obj.optnorm = 'l2';
                case 'linf'
                    obj.optnorm = 'linf';
                otherwise
                    error('Unknown atlas algorithm: set 1d or kd as atlas algorithm');
            end
        end        
        
        % convert fnl from tensor format to multiindex
        function fnlTensor2Multi(obj)
            fnl = obj.system.fnl;
            y = tensor_to_multi_index(fnl);
            obj.multiFnl = y;
        end

        
        % vector field compatiable to coco - autonomous case
        f = aut_ode(obj, x, p, data)
        
        % vector field compatible to coco - nonautonomous case
        f = ode_het(obj, t, x, p, data);
        
        % Jacobian of ode_het w.r.t x
        J = ode_het_dx(obj, t, x, p, data);

        % Jacobian of ode_het w.r.t t
        J = ode_het_dp(obj, t, x, p, data);

        % Jacobian of ode_het w.r.t t
        J = ode_het_dt(obj, t, x, p, data);
        
        % internal force in M\ddot{u}+N(u,\dot{u}) = F(t,p)
        y = Nhan(obj,u,v);
        
        % derivative of N w.r.t u
        y = dNdu(obj,u,v);
        
        % derivative of N w.r.t v
        y = dNdv(obj,u,v);
        
        % external force in M\ddot{u}+N(u,\dot{u}) = F(t,p)
        y = Fext(obj,t,p);
        
        % derivative of Fext w.r.t p
        y = dFextdp(obj,t,p);
        
        % setup coco options
        prob = cocoSet(obj, prob)
        
        % extract backbone curve for mode whose natural frequency is
        % closest to omega among all natural frequencies
        bd = extract_backbone(obj, omega, varargin) 
        
        % extract Force Response Curve for given frequency range
        bd = extract_FRC(obj, parName, parRange, varargin)
        
        % extract Force Response Curve for given frequency range starting
        % from an initial periodic orbit
        bd = extract_FRC_fromIC(obj, omega_range,IC,omegaIC ,varargin)

        % extract FRC using forward simulation
        bd = forward_FRC(obj, omega_range, varargin)
        
        % continuation of TR bifurcation periodic orbit
        bd = FRC_po2TR(obj, runid, lab, parRange, varargin)
        
        % continuation of tori born from a TR bifurcation periodic orbit
        bd = FRC_TR2tor(obj, runid, lab, omega_range, nseg, epsilon, varargin)
        
        % continuation of tori from a initial solution guess
        bd = FRC_isol2tor(obj, omega_range, t0, x0, p0);
        
        % extract Stability Diagram from full system
        SD = extract_Stability_Diagram(obj,omRange,epsRange,parName,p0,varargin)

        % extraction of damped backbone curve - optimization based
        % computation
        bds = extract_damped_backbone(obj,parRange,optdof,varargin)
                              
        % Coco PO Sweeps
        bds = coco_poSweeps(obj,oid,epSamp,omRange,varargin)
    end
end

