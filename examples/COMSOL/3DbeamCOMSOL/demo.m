%% Finding a 2D SSM for a 3D finite element beam
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using all phase space variables measurements. In this example, we consider the 
% clamped-clamped 3D beam that is modelled in COMSOL Multiphysics.
clear all
close all

% First we need to establish the connection between COMSOL and MATLAB from
% COMSOL server. See details in: https://doc.comsol.com/5.6/doc/com.comsol.help.llmatlab/llmatlab_ug_start.5.04.html
run ../../../install.m
%% generate model with a given point load
f0 = 1e5;

[model, M, C, K, f_0, Null, Nullf, ud, outdof, u0, disp_idx] = build_model(f0);
%%
% Construction of M,C,K of the finite element model in COMSOL. We use the
% function provided in COMSOL to extract the internal stiffness force.
[fint,dfint] = get_fint(K,model,Null, Nullf, ud, u0, disp_idx);

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
set(DS.Options,'outDOF',outdof);

%%
% We assume periodic forcing of the form
% 
% $$\mathbf{f}^{ext}(\phi) = \mathbf{f}_0\cos(\phi)=\frac{\mathbf{f}_0}{2}e^{i\phi} 
% + \frac{\mathbf{f}_0}{2}e^{-i\phi}  $$
% 
% Fourier coefficients of Forcing

epsilon = 0.01;
kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
%% Linear Modal analysis and SSM setup
[V,D,W] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 1,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',false,'COMPtype','second')
set(S.FRCOptions, 'nt', 2^7, 'nRho', 300, 'nPar', 150, 'nPsi', 200, 'rhoScale', 5 )
set(S.FRCOptions, 'method','continuation ep' ) 
set(S.FRCOptions, 'outdof',outdof)

masterModes = [1,2]; 
S.choose_E(masterModes);
%% 
% choose frequency range around the first natural frequency
omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.95 1.05];
%% 
% extract forced response curve
order = 3;
FRC = S.extract_FRC('freq',omegaRange,order);
