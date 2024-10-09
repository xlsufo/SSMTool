function [R1,W1] = nonAut_assembleCoefficients(W1,W1_ik,R1,R1_ik,i,k,l)
% NONAUT_ASSEMBLECOEFFICIENTS This function inserts the computed

% non-autonomous SSM and RD coefficients into the respective datastructs
%
% [R1,W1] = NONAUT_ASSEMBLECOEFFICIENTS(W1,W1_ik,R1,R1_ik,i,k,l)
%
% W1:       non-autonomous SSM coefficient struct
% W1_ik:    non-autonomous SSM coefficients at order k for harmonic i
% R1:       non-autonomous RD coefficient struct
% R1_ik:    non-autonomous RD coefficients at order k for harmonic i
% i:        index of harmonic for which cohomol. eq. is solved
% k:        current order of SSM computation
% l:        dimension of SSM
%
% W1:       non-autonomous SSM coefficient struct
% R1:       non-autonomous RD coefficient struct
%
% See also: NONAUT_1STORDER_HIGHTERMS, NONAUT_2NDORDER_HIGHTERMS

R1(i).R(k+1).coeffs = R1_ik;
W1(i).W(k+1).coeffs = W1_ik;

if l > 1
    R1(i).R(k+1).ind    = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
    W1(i).W(k+1).ind    = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
    
else
    R1(i).R(k+1).ind = k;
    W1(i).W(k+1).ind = k;
end
end