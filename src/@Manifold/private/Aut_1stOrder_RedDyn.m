function [R_0i,RHS] = Aut_1stOrder_RedDyn(z_k, Lambda_M, Lambda_K, W_M, reltol,RHS,V_M,B)
% AUT_2NDORDER_REDDYN
%
% This function computes the kernel of the coefficient matrix for eigenvalue 
% pairs that are in resonance as described in the document ''Analytic Kernel of 
% Coefficient matrix and Symmetries using multi-index notation''.
% It then computes the autonomous reduced dynamics coefficients at
% order k using the first order systems routine. 
%
% [Rk] = AUT_2NDORDER_REDDYN(I,F,THETA,PHI,C,Lambda_K,Lambda, M,Vm,Ym,l,z_k)
%
% z_k:      number of distinct multi-indices at order k
% Lambda_M: Vector containing master mode eigenvalues
% Lambda_K:  direct product of multi-indices at order K with master mode
%           vector
% reltol:   relative tolerance for detecting internal resonances
% RHS:      RHS of invariance equation
% V_M:      right eigenvectors of system
% B:        system matrix
%
% R0_i:     order i autonomous reduced dynamics
% RHS:      updated RHS of invariance equation
%
% See also: AUT_1STORDER_SSM

%% COEFF_MATR_KERNEL Explicit kernel-construction of the coefficient-matrix

%SSM dimension
l         = size(Lambda_M,1);
%Compare for all combinations if singularity occurs
Lambda_Ci = Lambda_M - Lambda_K; % column vector - row vector
%threshold below which resonance occurs
ref       = min(abs(Lambda_M));
abstol = reltol*ref;
%find eigenvalues that trigger resonance
[I,F]  = find(abs(Lambda_Ci) < abstol); % I for eigenvalue and F for combination
r_k = length(I);
if r_k
    innerresonance = 1;
    
    % create E_F, E_I
    E_F = sparse( F, (1:r_k).', true(r_k,1), z_k, r_k);
    E_I = sparse( I, (1:r_k).', true(r_k,1), l, r_k);
    
    % create K_k, G_k
    K_k = khatri_rao_product(E_F, W_M(:,I));
    G_k = khatri_rao_product(E_F, E_I)';
else
    innerresonance = 0;
    
    K_k=[];
    G_k=[];
end

if innerresonance
    
    Skron = kron(speye(z_k),B*V_M);
    R_0i   = G_k.' * K_k' * RHS;
    RHS   = RHS - Skron*R_0i;
else
    R_0i   = sparse(l*z_k,1);
end

end