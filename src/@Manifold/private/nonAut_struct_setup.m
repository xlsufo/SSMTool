function [W1,R1,kappas,Fext_ord] = nonAut_struct_setup(l,N,nKappa,order,FEXT)
% NONAUT_STRUCT_SETUP This function is used to initialise the data structures 
%
% that are used in the non-autonomous SSM computation.
%
% [W1,R1,kappas,Fext_ord]= NONAUT_STRUCT_SETUP(l,N,nKappa,order,FEXT)
%
% l:        dimension of SSM
% N:        dimension of full phase space
% nKappa:   number of harmonics in forcing
% order:    approximation order up until which SSM is computed
% FEXT:     struct array containing the external forcing
%
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
% kappas:   array containing the harmonics that are present in the forcing
% Fext_ord: nonlinear order of the external forces for each harmonic
%
% See also: NONAUT_2NDORDER_WHISKER, NONAUT_1ST_ORDER_WHISKER

k_kappa   = size(FEXT.data(1).kappa,1);
kappas    = zeros(k_kappa,nKappa);
% intitalise data structures to store coefficients
idle = repmat(struct('coeffs',[],'ind',[]),order+1  , 1);
W1  = repmat(struct('kappa' ,[],'W',idle),nKappa, 1);
R1  = repmat(struct('kappa' ,[],'R',idle),nKappa, 1);

Fext_ord = zeros(1,nKappa);

for i = 1:nKappa
    Fext_ord(i)  = numel(FEXT.data(i).F_n_k);
    kappa = FEXT.data(i).kappa;
    W1(i).kappa = kappa;
    R1(i).kappa = kappa;
    kappas(:,i)  = kappa;
    
    W1(i).W(1).coeffs = sparse(N,1);
    W1(i).W(1).ind    = sparse(l,1);
    R1(i).R(1).coeffs = sparse(l,1);
    R1(i).R(1).ind    = sparse(l,1);
end
end