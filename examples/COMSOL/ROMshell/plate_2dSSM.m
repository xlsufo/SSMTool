%% Finding a 2D SSM for a 3D finite element Plate
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using all phase space variables measurements. In this example, we consider the 
% perforated cover plate modelled in COMSOL Multiphysics. Details of this
% geometry can be found in:
% https://doi.org/10.1016/j.jsv.2016.12.037

clear all
close all
run ../../../install.m
%%
% First we need to establish the connection between COMSOL and MATLAB from
% COMSOL server. See details in: https://doc.comsol.com/5.6/doc/com.comsol.help.llmatlab/llmatlab_ug_start.5.04.html
%% generate model from COMSOL
[model, M, C, K, Null, Nullf, ud, outdof] = build_model();
%%
% Construction of M,C,K of the finite element model in COMSOL. We use the
% function provided in COMSOL to extract the internal stiffness force.
[fint,dfint] = get_fint(K,C,model,Null, Nullf, ud);

n = length(M);
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
set(DS,'M',M,'C',C,'K',K);
set(DS.Options, 'Intrusion', 'none')

set(DS,'fnl_non',fint,'dfnl_non',dfint);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
%%
% We assume periodic forcing of the form
% 
% $$\mathbf{f}^{ext}(\phi) = \mathbf{f}_0\cos(\phi)=\frac{\mathbf{f}_0}{2}e^{i\phi} 
% + \frac{\mathbf{f}_0}{2}e^{-i\phi}  $$
% 
% Fourier coefficients of Forcing
f_0 = zeros(n,1);
f_0(outdof) = 1;
epsilon = 1.6;                                                                                              
masterModes = [1,2];
kappas = [-1; 1];
coeffs = [f_0 f_0]/2;

DS.add_forcing(coeffs, kappas,epsilon);
%% Linear Modal analysis and SSM setup
[V,D,W] = DS.linear_spectral_analysis();
%% SSM setup
% *Choose Master subspace (perform resonance analysis)*
S = SSM(DS);
S.choose_E(masterModes);
% %%
% setup options
set(S.Options, 'reltol', 1,'IRtol',0.1,'notation', 'multiindex','contribNonAuto',false)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 300, 'nPar', 100, 'nPsi', 100, 'rhoScale', 2 )
set(S.FRCOptions, 'method', 'level set') % 
set(S.FRCOptions, 'outdof',outdof)
set(DS.Options, 'outDOF',outdof)

% %% 
% choose frequency range around the first natural frequency
omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.95 1.05];% 
order = 7;
%% compute SSM coefficients and reduced dynamics
[W_0, R_0] = S.compute_whisker(order);

save('SSMCoeff.mat','W_0','R_0')