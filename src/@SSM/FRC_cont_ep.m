function varargout = FRC_cont_ep(obj,i,resModes,ORDER,mFreqs,parName,parRange)
% FRC_CONT_EP This function performs continuation of equilibrium points of
% slow dynamics. Each equilibirum point corresponds to a periodic orbit in
% the regular time dynamics. The continuation here starts from the guess of
% initial solution.
%
% FRC = FRC_CONT_EP(OBJ,OID,RESMODES,ORDER,MFREQS,PARNAME,PARRANGE)
%
% i:        index of subinterval
% resModes: master subspace
% order:    expansion order of SSM
% mFreqs:   internal resonance relation vector
% parName:  amp/freq continuation parameter
% parRange: continuation domain of parameter, which should be near the
%           value of natural frequency with index 1 in the mFreq if the continuation
%           parameter is freq


m = numel(mFreqs);
assert(numel(resModes)==2*m, 'The master subspace is not %dD.',2*m);

% get options
[nt,nPar,nCycle,sampStyle,initialSolver,outdof,saveIC,coordinates,p0,z0] = ...
    deal(obj.FRCOptions.nt, obj.FRCOptions.nPar, obj.FRCOptions.nCycle, ...
    obj.FRCOptions.sampStyle, obj.FRCOptions.initialSolver, ...
    obj.FRCOptions.outdof, obj.FRCOptions.saveIC, ...
    obj.FRCOptions.coordinates, obj.FRCOptions.p0, obj.FRCOptions.z0);

%% Checking whether internal resonance indeed happens
if isempty(obj.System.spectrum)
    [~,~,~] = obj.System.linear_spectral_analysis();
end
% eigenvalues Lambda is sorted in descending order of real parts
% positive imaginary part is placed first in each complex pair

lambda = obj.System.spectrum.Lambda(resModes);
lambdaRe = real(lambda);
lambdaIm = imag(lambda);
% check spectrum
check_spectrum(lambdaRe,lambdaIm,mFreqs);

%% SSM computation of autonomous part
obj.choose_E(resModes)
% compute autonomous SSM coefficients
ORDER = sort(ORDER);
[W_0,R_0] = obj.compute_whisker(ORDER(end));

for jj = 1:numel(ORDER)
order = ORDER(jj);
    
oid  = ['freqSubint',num2str(i),'Order',num2str(order)];
if jj>1 % take lowest order solution as initial guess
    old_oid  = ['freqSubint',num2str(i),'Order',num2str(ORDER(1))];

    sol_jminus1 = ep_read_solution('', [old_oid,'.ep'], 1);
    p0 = sol_jminus1.p;
    z0 = sol_jminus1.x;
end

% extract ind and coeff of expansion in reduced dynamics and check
% consistency. Here beta is coeff and kappa denotes ind
[beta,kappa] = extract_beta_kappa(m,order,R_0,mFreqs);

%% SSM computation of non-autonomous part
[iNonauto, rNonauto, kNonauto, W_1] = extract_red_nonauto(obj,mFreqs,lambdaIm,order,W_0,R_0,parName);
% Here iNonauto,rNonauto and kNonauto gives the indices for resonant happens
% value of leading order contribution and (pos) kappa indices with resonance

%% Construct COCO-compatible vector field
fdata = struct();
fdata.beta  = beta;
fdata.kappa = kappa;
fdata.lamdRe = lambdaRe(1:2:end-1);
fdata.lamdIm = lambdaIm(1:2:end-1);
fdata.mFreqs = mFreqs;
fdata.iNonauto = iNonauto;
fdata.rNonauto = rNonauto;
fdata.kNonauto = kNonauto;
% put W_0 and W_1 in fdata is a bad idea because it will be stored in disk
% for each saved continuation solution. As an alternative, we save W_0 and
% W_1 in disk here under folder data. When needed, they will be loaded into
% memory.
% fdata.W_0   = W_0;
% fdata.W_1   = W_1;
data_dir = fullfile(pwd,'data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
wdir = fullfile(data_dir,'SSM.mat');
SSMcoeffs = struct();
SSMcoeffs.W_0 = W_0(1:order);
SSMcoeffs.W_1 = W_1;
save(wdir, 'SSMcoeffs');
fdata.order = order;
fdata.modes = resModes;

ispolar = strcmp(coordinates, 'polar');
fdata.ispolar = ispolar;
fdata.isbaseForce = obj.System.Options.BaseExcitation;
if ispolar
    odefun = @(z,p) ode_2mDSSM_polar(z,p,fdata);
else
    odefun = @(z,p) ode_2mDSSM_cartesian(z,p,fdata);
end
funcs  = {odefun};

%% continuation of reduced dynamics w.r.t. parName
prob = coco_prob();
prob = cocoSet(obj.contOptions, prob);
if isempty(p0)
    if ~isempty(obj.System.fext)
        p0 = [obj.System.Omega; obj.System.fext.epsilon];
    else
        p0 = [obj.System.Omega; obj.System.Fext.epsilon];
    end
end
if isempty(z0)
    if ispolar
        z0 = 0.1*ones(2*m,1);
    else
        z0 = zeros(2*m,1);
    end
end

% construct initial guess equilibrium points
z0 = get_initial_sol(z0,p0,initialSolver,odefun,nCycle,ispolar);

% call ep-toolbox
prob = ode_isol2ep(prob, '', funcs{:}, z0, {'om' 'eps'}, p0);

% define monitor functions to state variables - depends on ispolar, the
% output args1/2 could be rho/theta or Real/Imag(z)
[prob,args1,args2] = monitor_states(prob,ispolar,m);

switch parName
    case 'freq'
        isomega = true;
        cont_args = [{'om'},args1(:)' ,args2(:)',{'eps'}];
    case 'amp'
        isomega = false;
        cont_args = [{'eps'},args1(:)' ,args2(:)',{'om'}];
    otherwise
        error('Continuation parameter should be freq or amp');
end
% add events in the case of uniform sampling
if strcmp(sampStyle, 'uniform')
    if isomega
        omSamp = linspace(parRange(1),parRange(2), nPar);
        prob   = coco_add_event(prob, 'UZ', 'om', omSamp);
    else
        epSamp = linspace(parRange(1),parRange(2), nPar);
        prob   = coco_add_event(prob, 'UZ', 'eps', epSamp);
    end
end

runid = coco_get_id(oid, 'ep');

fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  runid);

coco(prob, runid, [], 1, cont_args, parRange);

%% extract results of reduced dynamics at sampled frequencies
FRC = ep_reduced_results(runid,sampStyle,ispolar,isomega,args1,args2,'plot-off');

%% FRC in physical domain
om   = FRC.om;
epsf = FRC.ep;
state   = FRC.z;
numOm   = numel(om);
noutdof = numel(outdof);
Zout_frc  = cell(numOm,1);
Znorm_frc = zeros(numOm,1);
Aout_frc  = zeros(numOm,noutdof);

if saveIC
    Z0_frc = zeros(numOm,obj.System.N); % initial state
end
timeFRCPhysicsDomain = tic;

% Loop around a resonant mode
if obj.FRCOptions.nonAutoParRedCom
    nKappa  = obj.System.nKappa;
    kappa_set = zeros(nKappa,1);
    F_kappa = zeros(obj.dimSystem,nKappa); % each column corresponds to one kappa
    for ka_idx = 1:nKappa
        kappa_set(ka_idx) = obj.System.Fext.data(ka_idx).kappa;
        F_kappa(:,ka_idx) = obj.System.Fext.data(ka_idx).F_n_k(1).coeffs;
    end
    A = obj.System.A;           % A matrix
    B = obj.System.B;           % B matrix
    W_M = obj.E.adjointBasis;   % Right eigenvectors of the modal subspace
    V_M = obj.E.basis;          % Left eigenvectors of the modal subspace
    reltol = obj.Options.reltol;
    % *Near external resonances and Leading order reduced dynamics*
    Lambda_M_vector = obj.E.spectrum;

    % transfer minimum data for reducing communication load
    W1 = leading_order_nonauto_SSM(A,B,W_M,V_M,Lambda_M_vector,kappa_set,F_kappa,reltol,kNonauto,iNonauto,rNonauto,om);

    % map it back to physical coordinates
    for j=1:numel(om)
        % Forced response in Physical Coordinates
        statej = state(j,:);
        % create compatible W1j
        idle = repmat(struct('coeffs',[],'ind',[]),1  , 1);
        W1j  = repmat(struct('kappa' ,[],'W',idle),nKappa, 1);
        assert(norm(kappa_set-W1{j}.kappas)<1e-6,'the orders of kappa are not consistent');
        for i = 1:nKappa
            W1j(i).kappa = kappa_set(i);
            W1j(i).W(1).coeffs = W1{j}.coeffs(:,i);
            W1j(i).W(1).ind    = sparse(obj.dimManifold,1);
        end
        if obj.System.Options.BaseExcitation
            epsf(j) = epsf(j)*(om(j))^2;
        end
        [Aout, Zout, z_norm, Zic] = compute_full_response_2mD_ReIm(W_0(1:order), W1j, statej, epsf(j), nt, mFreqs, outdof);

        % collect output in array
        Aout_frc(j,:) = Aout;
        Zout_frc{j}   = Zout;
        Znorm_frc(j)  = z_norm;

        if saveIC
            Z0_frc(j,:) = Zic; % initial state
        end

    end

else
    parfor j = 1:numel(om)
        % compute non-autonomous SSM coefficients
        if isomega
            set_omega(obj,om(j));
            if obj.Options.contribNonAuto
                [W_1j, R_1] = obj.compute_perturbed_whisker(0,[],[]);
                
                R_10 = R_1(kNonauto).R.coeffs;
                assert(~isempty(R_10), 'Near resonance does not occur, you may tune tol');
                f = R_10(2*iNonauto-1);

                assert(norm(f-rNonauto)<1e-3*norm(f), 'inner resonance assumption does not hold');

    %             fprintf('the forcing frequency %.4d is nearly resonant with the eigenvalue %.4d + i%.4d\n', om(j), real(lambda(1)),imag(lambda(1)))
            else
                W_1j = [];
            end
        else
            W_1j = W_1;
        end
        % Forced response in Physical Coordinates
        if obj.System.Options.BaseExcitation
            epsf(j) = epsf(j)*(om(j))^2;
        end
        statej = state(j,:);
        [Aout, Zout, z_norm, Zic] = compute_full_response_2mD_ReIm(W_0(1:order), W_1j, statej, epsf(j), nt, mFreqs, outdof);

        % collect output in array
        Aout_frc(j,:) = Aout;
        Zout_frc{j}   = Zout;
        Znorm_frc(j)  = z_norm;

        if saveIC
            Z0_frc(j,:) = Zic; % initial state
        end
    end
end

maptime = toc(timeFRCPhysicsDomain);
fprintf(['Time spent on mapping the FRC to physical coordinates maptime ' num2str(maptime) ' seconds. \n'])
%%
% Record output
FRC.Aout_frc  = Aout_frc;
FRC.Zout_frc  = Zout_frc;
FRC.Znorm_frc = Znorm_frc;

if saveIC
FRC.Z0_frc    = Z0_frc; % initial state
end
FRCinfo = struct();
FRCinfo.timeFRCPhysicsDomain = toc(timeFRCPhysicsDomain);
FRCinfo.SSMorder   = order;
FRCinfo.SSMnonAuto = obj.Options.contribNonAuto;
FRCinfo.SSMispolar = ispolar;

% convert results to cell array
FRC = array2structArray(FRC,saveIC);

% Plot Plot FRC in system coordinates
% if isomega
%     plot_FRC_full(FRC,outdof,order,'freq','lines');
% else
%     plot_FRC_full(FRC,outdof,order,'freq','lines');
% end


FRCs{jj} = FRC;
end

varargout{1} = FRCs;
fdir = fullfile(data_dir,runid,'SSMep.mat');
save(fdir, 'FRC','FRCinfo');
end


function check_spectrum(lambdaRe,lambdaIm,mFreqs)

flags1 = abs(lambdaRe(1:2:end-1)-lambdaRe(2:2:end))<1e-6*abs(lambdaRe(1:2:end-1)); % same real parts
flags1 = all(flags1);
flags2 = abs(lambdaIm(1:2:end-1)+lambdaIm(2:2:end))<1e-6*abs(lambdaIm(1:2:end-1)); % opposite imag parts
flags2 = all(flags2);
freqs  = lambdaIm(1:2:end-1);
freqso = freqs - dot(freqs,mFreqs(:))*mFreqs(:)/sum(mFreqs.^2);
flags3 = norm(freqso)<0.1*norm(freqs);

assert(flags1, 'Real parts do not follow complex conjugate relation');
assert(flags2, 'Imaginary parts do not follow complex conjugate relation');
assert(flags3, 'Internal resonnace is not detected for given master subspace');

end


function [beta,kappa] = extract_beta_kappa(m,order,R_0,mFreqs)
beta  = cell(m,1); % coefficients - each cell corresponds to one mode
kappa = cell(m,1); % exponants
Em = eye(m);
for k = 2:order
    R = R_0(k);
    coeffs = R.coeffs;
    ind = R.ind;
    if ~isempty(coeffs)
        for i=1:m
            betai = coeffs(2*i-1,:);
            [~,ki,betai] = find(betai);
            kappai = ind(ki,:);
            % check resonant condition
            l = kappai(:,1:2:end-1);
            j = kappai(:,2:2:end);
            nk = numel(ki);
            rm = repmat(mFreqs(:)',[nk,1]);
            flagi = dot(l-j-repmat(Em(i,:),[nk,1]),rm,2);
            assert(all(flagi==0), 'Reduced dynamics is not consisent with desired IRs');
            % assemble terms
            beta{i}  = [beta{i} betai];
            kappa{i} = [kappa{i}; kappai];
        end
    end
end
end

function [iNonauto, rNonauto, kNonauto, W_1] = extract_red_nonauto(obj,mFreqs,lambdaIm,order,W0,R0,parName)
iNonauto = []; % indices for resonant happens
rNonauto = []; % value of leading order contribution
kNonauto = []; % (pos) kappa indices with resonance
kappa_set= obj.System.kappas; % each row corresponds to one kappa
kappa_pos = kappa_set(kappa_set>0);
num_kappa = numel(kappa_pos); % number of kappa pairs
for k=1:num_kappa
    kappak = kappa_pos(k);
    idm = find(mFreqs(:)==kappak); % idm could be vector if there are two frequencies are the same
    if strcmp(parName,'freq')
        obj.System.Omega = lambdaIm(2*idm(1)-1);
    end


    [W_1, R_1] = obj.compute_perturbed_whisker(0,[],[]);

    idk = find(kappa_set==kappak);
    R_10 = R_1(idk).R.coeffs;
    r = R_10(2*idm-1);

    iNonauto = [iNonauto; idm];
    rNonauto = [rNonauto; r];
    kNonauto = [kNonauto; idk];
end
end

function z0 = get_initial_sol(z0,p0,initialSolver,odefun,nCycle,ispolar)
switch initialSolver
    case 'fsolve'
        % fsolve to approximate equilibrium
        fsolveOptions = optimoptions('fsolve','MaxFunctionEvaluations',100000,...
            'MaxIterations',1000000,'FunctionTolerance', 1e-10,...
            'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-10);
        z0 = fsolve(@(z) odefun(z,p0),z0,fsolveOptions);
    case 'forward'
        % forward simulation to approach equilibirum
        tspan = [0 nCycle*2*pi/p0(1)]; %nCycle
        odefw = @(t,z,p) odefun(z,p);
        opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
        [~,y0] = ode45(@(t,y) odefw(t,y,p0), tspan, z0, opts);
        [~,y] = ode45(@(t,y) odefw(t,y,p0), [0 2*pi/p0(1)], y0(end,:));
        [~, warnId] = lastwarn;

        if any(isnan(y(:))) || strcmp(warnId,'MATLAB:ode45:IntegrationTolNotMet')
            warning('Reduced dynamics with IRs in polar form diverges with [0.1 0.1 0.1 0.1]');
        else
            z0 = y(end,:)';
        end
end

if ispolar % regularize initial solution if it is in polar form
    z0(2:2:end) = mod(z0(2:2:end),2*pi); % phase angles in [0,2pi]
    m = numel(z0)/2;
    for k=1:m
        if z0(2*k-1)<0
            z0(2*k-1) = -z0(2*k-1);      % positive amplitudes
            z0(2*k) = z0(2*k)+pi;
        end
    end
end
end

function [prob,varargout] = monitor_states(prob,ispolar,m)
if ispolar
    rhoargs = cell(m,1);
    thargs  = cell(m,1);
    for k=1:m
        rhoargs{k} = strcat('rho',num2str(k));
        thargs{k}  = strcat('th',num2str(k));
    end
    prob = coco_add_pars(prob, 'radius', 1:2:2*m-1, rhoargs(:)');
    prob = coco_add_pars(prob, 'angle', 2:2:2*m, thargs(:)');
    varargout{1} = rhoargs;
    varargout{2} = thargs;
else
    Reargs = cell(m,1);
    Imargs = cell(m,1);
    for k=1:m
        Reargs{k} = strcat('Rez',num2str(k));
        Imargs{k} = strcat('Imz',num2str(k));
    end
    prob = coco_add_pars(prob, 'realParts', 1:2:2*m-1, Reargs(:)');
    prob = coco_add_pars(prob, 'imagParts', 2:2:2*m, Imargs(:)');
    varargout{1} = Reargs;
    varargout{2} = Imargs;
end
end

function FRCy = array2structArray(FRCx,saveIC)

om    = FRCx.om;
state = FRCx.z;
rho   = FRCx.rho;
th    = FRCx.th;
ep    = FRCx.ep;
st    = FRCx.st;
SNidx = FRCx.SNidx;
HBidx = FRCx.HBidx;
Aout  = FRCx.Aout_frc;
Zout  = FRCx.Zout_frc;
Znorm = FRCx.Znorm_frc;
if saveIC
Zic   = FRCx.Z0_frc;
end

numPts = numel(om);
FRCy   = cell(numPts,1);
for i=1:numPts
    FRC = struct('rho', rho(i,:), 'stability', st(i), 'Omega', om(i) ,...
    'epsilon', ep(i), 'Aout', Aout(i,:), 'Zout', Zout{i},...
    'Znorm', Znorm(i), 'Zic', [], 'isSN', false, 'isHB', false, ...
    'state', state(i,:), 'th', th(i,:));
    if saveIC
        FRC.Zic = Zic(i,:);
    end
    if ismember(i,SNidx)
        FRC.isSN = true;
    end
    if ismember(i,HBidx)
        FRC.isHB = true;
    end
    FRCy{i} = FRC;
end

FRCy = cat(1,FRCy{:});
end 


function set_omega(obj,omega)
obj.System.Omega = omega;
end