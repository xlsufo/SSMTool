clc; addpath('utils') 
%% Load data
R_0 = load('SSMCoeff.mat','R_0');
R_0 = R_0.R_0;
W_0 = load('SSMCoeff.mat','W_0'); 
W_0 = W_0.W_0;
mFreqs = 1; order = 7;
% check reduced dynamics (consistent with expected internal resonance)
[beta,kappa] = check_auto_reduced_dynamics(R_0(1:order),order,mFreqs);
%% create DS and S

DSorder = 2;
DS = DynamicalSystem(DSorder);
set(DS,'M',M,'C',C,'K',K);
set(DS.Options, 'Intrusion', 'none')
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
[V,D,W] = DS.linear_spectral_analysis();
% external forcing
kappas  = [-1; 1];
coeffs  = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
% %% Linear Modal analysis and SSM setup
S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex','contribNonAuto',false);
resonant_modes = [1 2];
S.choose_E(resonant_modes);

%% SSM computation of non-autonomous part
lambda = S.System.spectrum.Lambda(resonant_modes);
lambdaRe = real(lambda);
lambdaIm = imag(lambda);
iNonauto = []; % indices for resonant happens
rNonauto = []; % value of leading order contribution
kNonauto = []; % (pos) kappa indices with resonance
kappa_set= [S.System.Fext.data.kappa]; % each row corresponds to one kappa
kappa_pos = kappa_set(kappa_set>0);
num_kappa = numel(kappa_pos); % number of kappa pairs
for k=1:num_kappa
    kappak = kappa_pos(k);
    idm = find(mFreqs(:)==kappak); % idm could be vector if there are two frequencies are the same
    S.System.Omega = lambdaIm(2*idm(1)-1);

    [W_1, R_1] = S.compute_perturbed_whisker(0,[],[]);

    idk = find(kappa_set==kappak);
    R_10 = R_1(idk).R.coeffs;
    r = R_10(2*idm-1);

    iNonauto = [iNonauto; idm];
    rNonauto = [rNonauto; r];
    kNonauto = [kNonauto; idk];
end
%% Construct COCO-compatible vector field
% create data to vector field
lamd  = struct();
lamd.lambdaRe = lambdaRe; lamd.lambdaIm = lambdaIm;
Nonauto = struct();
Nonauto.iNonauto = iNonauto; Nonauto.rNonauto = rNonauto; Nonauto.kNonauto = kNonauto;
[fdata,data_dir] = create_reduced_dynamics_data(beta,kappa,lamd,mFreqs,Nonauto,W_0(1:order),W_1,order,resonant_modes);

ispolar = strcmp(S.FRCOptions.coordinates, 'polar');
fdata.ispolar = ispolar;
fdata.isbaseForce = S.System.Options.BaseExcitation;

if ispolar
    odefun = @(z,p) ode_2mDSSM_polar(z,p,fdata);
else
    odefun = @(z,p) ode_2mDSSM_cartesian(z,p,fdata);
end
funcs  = {odefun};

%% continuation problem
prob = coco_prob();
set(S.contOptions,'PtMX',300,'h_max',2,'h_min',1e-3);
prob = cocoSet(S.contOptions, prob);
% construct initial solution
p0 = [S.System.Omega; S.System.fext.epsilon];
m  = numel(mFreqs);
nCycle = S.FRCOptions.nCycle;
[p0,z0] = initial_fixed_point(p0,S.FRCOptions.initialSolver,ispolar,...
    odefun,nCycle,m);

% call ep-toolbox
prob = ode_isol2ep(prob, '', funcs{:}, z0, {'om' 'eps'}, p0);
% define monitor functions to state variables
[prob, args1, args2] = monitor_states(prob, ispolar, m);

isomega = true;
runid = coco_get_id(['contep',num2str(order)], 'ep');
fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  runid);

cont_args = [{'om'},args1(:)' ,args2(:)',{'eps'}];

omega0 = imag(S.E.spectrum(1));
omegaRange = [1.208413758931806e+03 1.335615207240417e+03];% 
coco(prob, runid, [], 1, cont_args, omegaRange);

%% extract results of reduced dynamics at sampled frequencies
FRC = ep_reduced_results(runid,S.FRCOptions.sampStyle,ispolar,isomega,args1,args2);

%% FRC in physical domain
ns = numel(FRC.om);
% time-dependent leading-order contribution
if S.Options.contribNonAuto
    W_1 = cell(ns,1);
    for j=1:ns
        j/ns
        S.System.Omega = FRC.om(j);
        if isomega
            [W1, R1] = S.compute_perturbed_whisker(0,[],[]);
    
            R10 = R1(kNonauto).R.coeffs;
            assert(~isempty(R10), 'Near resonance does not occur, you may tune tol');
            f = R10(2*iNonauto-1);
    
            assert(norm(f-rNonauto)<1e-3*norm(f), 'inner resonance assumption does not hold');
        end
        W_1{j} = W1;
    end
else
    W_1 = [];
end
FRCdata = struct();        FRCdata.isomega = isomega;
FRCdata.mFreqs  = mFreqs;  FRCdata.order  = order;
FRCdata.ispolar = ispolar; FRCdata.modes  = resonant_modes;
outdof = 2809;
FRC = FRC_reduced_to_full(S.FRCOptions.nt,Nonauto,S.Options.contribNonAuto,FRC,FRCdata,W_0,W_1,outdof);
%% Plot Plot FRC in system coordinates
plot_frc_full(FRC.om,FRC.Znorm_frc,outdof,FRC.Aout_frc,FRC.st,order,'freq','lines',{FRC.SNidx,FRC.HBidx});

%% backbone curve
gamma = compute_gamma(R_0);
order = [3,5,7];
set(S.FRCOptions,'rhoScale',1.5);
f1 = figure('Name','Norm');
if isnumeric(outdof)
    f2 = figure('Name',['Amplitude at DOFs ' num2str(outdof(:)')]);
else
    f2 = figure('Name','Amplitude at DOFs');
end
figs = [f1, f2];
colors = get(0,'defaultaxescolororder');
f3 = figure('Name','Amplitude in rho');

for k=1:numel(order)
    % compute backbone
    gidx  = round((order(k)+0.5)/2)-1;
    rho = compute_rho_grid(omegaRange,S.FRCOptions.nPar,...
        S.FRCOptions.rhoScale,gamma(1:gidx),lambda(1),S.FRCOptions.nRho);
    
    [~,b] = frc_ab(rho, 0, gamma(1:gidx), lambda(1));
    omega = b./rho;
    idx = [find(omega<omegaRange(1)) find(omega>omegaRange(2))];
    rho(idx) = []; omega(idx) = [];
    
    % Backbone curves in Physical Coordinates
    stability = true(size(rho)); psi = zeros(size(rho)); epsilon = 0;
    BC = compute_output_polar2D(rho,psi,stability,epsilon,omega,W_0(1:order(k)),[],1,S.FRCOptions.nt, false, outdof);
    % plotting
    nomegas = numel([BC.Omega]);
    for idx = 1:nomegas
        BC(idx).Omega = BC(idx).Omega;
    end
    plot_FRC(BC,outdof,order(k),'freq','lines',figs,colors(k,:));
%     xlim([1.5e5 1.7e5]/2/pi); ylim([0 4e-6]);
    figure(f3); hold on
    plot(omega/(2*pi),rho,'Color',colors(k,:));
end
