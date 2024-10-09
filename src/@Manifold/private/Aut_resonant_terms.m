function [I,F] = Aut_resonant_terms( Lambda_M, Lambda_K, reltol)
% AUT_RESONANT_TERMS  This function finds the combinations of multi-indices
%
% at order k and master mode eigenvalues that lead to internal resonances.
%
% [I,F] = AUT_RESONANT_TERMS( Lambda_M_vector, Lambda_Mk_vector, reltol)
%
% Lambda_M: vector containing master mode eigenvalues          
% Lambda_K: direct product of multi-indices at order k with master mode
%           vector           
% reltol:   relative tolerance for internal resonance
%
% I:        indices of the master modes that lead to resonances        
% F:        indeices of multi-indices that lead to resonance
%
% See also: NONAUT_RESONANT_TERMS

%Compare for all combinations if singularity occurs
Lambda_Ci = Lambda_M - Lambda_K; % column vector - row vector
%threshold below which resonance occurs
ref       = min(abs(Lambda_M));
abstol = reltol*ref;
%find eigenvalues that trigger resonance
[I,F]  = find(abs(Lambda_Ci) < abstol); % I for eigenvalue and F for combination
end