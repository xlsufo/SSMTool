%% Shallow-curved shell structure with geometric nonlinearities
% Finite element model used in the following reference:
% 
% Jain, S., & Tiso, P. (2018). Simulation-free hyper-reduction for geometrically 
% nonlinear structural dynamics: a quadratic manifold lifting approach. _Journal 
% of Computational and Nonlinear Dynamics_, _13_(7), 071003. <https://doi.org/10.1115/1.4040021 
% https://doi.org/10.1115/1.4040021>
% 
% Finite element code taken from the following package:
% 
% Jain, S., Marconi, J., Tiso P. (2020). YetAnotherFEcode (Version v1.1). Zenodo. 
% <http://doi.org/10.5281/zenodo.4011282 http://doi.org/10.5281/zenodo.4011282>
% 
% 
%% 
% *System parameters*

clear all
nDiscretization = 10; % Discretization parameter (#DOFs is proportional to the square of this number)
epsilon = 0.1; % converge at order 5
%% generate model

[M,C,K,fnl,f_0,outdof] = build_model(nDiscretization);
n = length(M); % number of degrees of freedom
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
%% Dynamical system setup 
% We consider the forced system
% 
% $$\mathbf{M}\ddot{\mathbf{x}}+\mathbf{C}\dot{\mathbf{x}}+\mathbf{K}\mathbf{x}+\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})=\epsilon\mathbf{f}^{ext}(\mathbf{\Omega}t),$$
% 
% which can be written in the first-order form as 
% 
% $$\mathbf{B}\dot{\mathbf{z}}	=\mathbf{Az}+\mathbf{F}(\mathbf{z})+\epsilon\mathbf{F}^{ext}(\mathbf{\phi}),\\\dot{\mathbf{\phi}}	
% =\mathbf{\Omega}$$
% 
% where
% 
% $\mathbf{z}=\left[\begin{array}{c}\mathbf{x}\\\dot{\mathbf{x}}\end{array}\right],\quad\mathbf{A}=\left[\begin{array}{cc}-\mathbf{K} 
% & \mathbf{0}\\\mathbf{0} & \mathbf{M}\end{array}\right],\mathbf{B}=\left[\begin{array}{cc}\mathbf{C} 
% & \mathbf{M}\\\mathbf{M} & \mathbf{0}\end{array}\right],\quad\quad\mathbf{F}(\mathbf{z})=\left[\begin{array}{c}\mathbf{-\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})}\\\mathbf{0}\end{array}\right],\quad\mathbf{F}^{ext}(\mathbf{z},\mathbf{\phi})=\left[\begin{array}{c}\mathbf{f}^{ext}(\mathbf{\phi})\\\mathbf{0}\end{array}\right]$.

DSorder = 2;
DS = DynamicalSystem(DSorder);
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
%% 
% We assume periodic forcing of the form
% 
% $$\mathbf{f}^{ext}(\phi) = \mathbf{f}_0\cos(\phi)=\frac{\mathbf{f}_0}{2}e^{i\phi} 
% + \frac{\mathbf{f}_0}{2}e^{-i\phi}  $$
% 
% Fourier coefficients of Forcing

kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
%% Linear Modal analysis and SSM setup

[V,D,W] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2,3,4]; 
S.choose_E(masterModes);
%% Forced response curves using SSMs
% Obtaining *forced response curve* in reduced-polar coordinate

order = 3;
set(S.Options, 'reltol', 0.5,'IRtol',0.05,'notation', 'multiindex','contribNonAuto',true)
%% 
% choose frequency range around the first natural frequency

set(S.FRCOptions,'coordinates','polar','initialSolver','forward');
set(S.contOptions, 'h_min', 1e-2,'h_max',2,'PtMX',300);
omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.92 1.07];
mFreq = [1 2];
%% 
% extract forced response curve

p0 = [omegaRange(1) epsilon]';
z0 = 1e-3*[1 1 1 1]';
S.SSM_isol2ep('isol-3',masterModes,order,mFreq,'freq',omegaRange,outdof,{p0,z0});
%%
% increase order
order = 5;
sol = ep_read_solution('','isol-3.ep',1);
set(S.FRCOptions,'initialSolver','fsolve');
S.SSM_isol2ep('isol-5',masterModes,order,mFreq,'freq',omegaRange,outdof,{sol.p,sol.x});