function bd = forward_FRC(obj, omega_range, varargin)
% FORWARD_FRC This function extract FRC of a second-order system using
% shooting method combined with parameter continuation.
% 
% BD = FORWARD_FRC(OBJ, OMEGA_RANGE, VARARGIN)
%
% OBJ:         A coco object with dynamical system included
% OMEGA_RANGE: Continuation domain of excitation frequency
% VARARGIN:    {omega0, x0} - initial solutoin
%
% See also: EXTRACT_FRC


%% initial setup
dir_name = obj.Options.dir_name;
N = obj.system.N;
n = obj.system.n;
if isempty(varargin)
    omega0 = omega_range(1); 
else
    omega0 = varargin{1}{1};
end
T       = 2*pi/omega0*obj.periodsRatio;
tf      = obj.nCycles*T;
outdof  = obj.outdof;
epsilon = obj.system.fext.epsilon;
p0 = [omega0;epsilon];

if ~isempty(obj.system.Omega)
    assert(numel(obj.system.Omega)==1, 'coco run assumes single freq component');
end
assert(obj.system.order == 2, 'fnl avaliable only for second-order systems')

opts = struct();
opts.ItMX   = obj.Options.IntItMX;
opts.RelTol = obj.Options.RelTol;
opts.Nsteps = obj.Options.Nsteps;
opts.alpha  = obj.Options.alpha;

if isempty(varargin)
    switch obj.initialGuess
        case 'forward'
            %% initial solution by forward simulation
            % ode45 is used here. Integration option may be added in future
            x0 = zeros(N,1);
            opts.Nsteps = obj.nCycles*obj.Options.Nsteps;
            [uend,vend] = Newmark(obj.system.M,@obj.Nhan,@obj.dNdu,...
                @obj.dNdv,@obj.Fext,@obj.dFextdp,0,tf,x0,p0,opts);         % transient
            x0 = [uend;vend];
            opts.Nsteps = obj.Options.Nsteps;
            [uend,vend] = Newmark(obj.system.M,@obj.Nhan,@obj.dNdu,...
                @obj.dNdv,@obj.Fext,@obj.dFextdp,0,T,x0,p0,opts);          % steady
            x1 = [uend;vend];
        case 'linear'
            %% initial solution by solving linear equation M\ddot{x}+C\dot{x}+Kx = F cos(O*t)
            mass = obj.system.M;
            damp = obj.system.C;
            stif = obj.system.K;
            fext = obj.system.fext.data.f_n_k;
            fext = 2*fext.coeffs;
            fext = epsilon*fext;
            kapa = obj.system.fext.data.kappa; kapa = kapa(1);
            if kapa<0
                kapa = -kapa;
            end
            qcom = (-(kapa*omega0)^2*mass+1i*kapa*omega0*damp+stif)\fext(:,1);
            t0   = [0,T];
            solx = qcom*exp(1i * kapa * omega0 * t0);
            solv = 1i * kapa * omega0 *solx;
            xx   = real([solx; solv]);  
            x0   = xx(:,1);
            x1   = xx(:,2);
    end
else
    x0 = varargin{1}{2};
    x1 = x0;
end

%% continuation excitation frequency
% setup coco
prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'forward', 'autonomous', false, 'order', 2);
switch lower(obj.forwardSolver)
    case 'newmark'
        prob = coco_set(prob, 'forward', 'ODEsolver', @Newmark);
    case 'galpha'
        prob = coco_set(prob, 'forward', 'ODEsolver', @Galpha);
end
forward_args = {obj.system.M, @obj.Nhan, @obj.dNdu, @obj.dNdv, @obj.Fext,...
    [], @obj.dFextdp, 0, T, x0, x1, {'omega' 'eps'}, p0};
prob = ode_isol2forward(prob, '',forward_args{:});

[data, uidx] = coco_get_func_data(prob, 'forward', 'data', 'uidx');
poData = struct();
poData.periodsRatio = obj.periodsRatio;
poData.dim = N;
bc_funcs = {@po_bc, @po_bc_du};
prob = coco_add_func(prob, 'bc', bc_funcs{:}, poData, 'zero', 'uidx', ...
  uidx([data.x0_idx, data.x1_idx, data.T0_idx, data.T_idx, data.p_idx(1)]));

% track amplitude of outdof AND stability
ampData = struct();
ampData.dof  = outdof;
ampData.dim  = N;
ampData.M    = obj.system.M;
ampData.Nhan = @obj.Nhan;
ampData.dNdu = @obj.dNdu;
ampData.dNdv = @obj.dNdv;
ampData.Fext = @obj.Fext;
ampData.dFextdp = @obj.dFextdp;
ampData.opts = opts;
numoutdof    = numel(outdof);
ampNames     = cell(1, numoutdof+1);
for k = 1:numoutdof
   ampNames{k} = strcat('amp',num2str(outdof(k))); 
end
ampNames{k+1} = 'stab';
prob = coco_add_func(prob, 'amp', @amplitude, ampData, 'regular', ampNames,...
    'uidx', uidx([data.x0_idx, data.T0_idx, data.T_idx, data.p_idx]));

cont_args = {1, [{'omega'} {'eps'} ampNames(:)'], [omega_range(1) omega_range(end)]};
    
runid = coco_get_id(dir_name, 'forwardFRC');
fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  runid);

bd = coco(prob, runid, [], cont_args{:});

end


function [data, y] = po_bc(~, data, u)

N  = data.dim;
x0 = u(1:N);
x1 = u(N+1:2*N);
T0 = u(2*N+1);
T  = u(2*N+2);
om = u(2*N+3);

y = [x1-x0; T0; T-2*pi*data.periodsRatio/om];

end

function [data, J] = po_bc_du(~, data, u)

N  = data.dim;
om = u(end);
J  = zeros(N+2,2*N+3);
J(1:N,1:N) = -eye(N);
J(1:N,N+1:2*N) = eye(N);
J(N+1,2*N+1) = 1;
J(N+2,2*N+2) = 1;
J(N+2,2*N+3) = 2*pi/om^2;

% [~,Jp] = coco_ezDFDX('f(o,d,x)',prob,data,@po_bc,u);

end

function [data, y] = amplitude(~, data, u)

N  = data.dim;
x0 = u(1:N);
T0 = u(N+1);
T  = u(N+2);
p  = u(N+3:end);
% simply evaluate amplitudes at outdofs
% [~,~,zt] = Newmark(data.M,data.Nhan,data.dNdu,...
%     data.dNdv,data.Fext,data.dFextdp,T0,T,x0,p,data.opts,data.dof);
% y = max(abs(zt),[],2);

% evaluate amplitudes at outdofs AND stability of periodic orbit
[~,~,Jup,Jvp,zt] = Newmark(data.M, data.Nhan, data.dNdu, data.dNdv, ...
    data.Fext, data.dFextdp, T0, T, x0, p, data.opts, 'var', data.dof);
% rewrite results in first-order form
y = max(abs(zt),[],2);
Jx0     = [Jup(:,2:data.dim+1); Jvp(:,2:data.dim+1)];
leadEig = eigs(Jx0,1);
stab    = double(abs(leadEig)<1);

y = [y; stab];
end


