function [WR]                    = nonAut_W1R0_plus_W0R1(data,k,W0,W1,R0,R1)
% NONAUT_W1R0_PLUS_W0R1 This function computes the multiplication of SSM 
%
% coefficients and reduced dynamics coefficients, collecting all resulting 
% terms that correspond to the set of multi-indices of order k. These
% compositions are first order terms, so it is either non-aut Red. Dyn and
% aut. SSM coefficients or vice versa that are multiplied.
%
% [Fk]= NONAUT_W1R0_PLUS_W0R1(F,W,nl_order,K,data)
%
% R0:       autonomous Reduced Dynamics Coefficients
% W0:       autonomous SSM coefficients
% R1:       non-autonomous Reduced Dynamics Coefficients
% W1:       non-autonomous SSM coefficients
% k:        order of computation
% data:     struct containing necessary information for the computation
%
% WR:       Composed terms W1R0 + W0R1 |_k
%
% See also: COEFFS_MIXED_TERMS, NONAUT_2NDORDER_HIGHTERMS, NONAUT_1STORDER_HIGHTERMS

z_k = nchoosek(k+data.l-1,data.l-1);

W1R0 = sparse(data.N,z_k);
W0R1 = sparse(data.N,z_k);

% Terms with non_aut SSM coefficients
data.mix = 'W1';
for m = 1:k % includes the zeroth order of W1, no order k non-aut SSM coeffs W1(k+1)
    if  ~isempty(W1.W(m).coeffs)
        W1R0 = W1R0 + coeffs_mixed_terms(k,m, W1.W, R0,data,'W1');
    end
end

% Terms with non_aut reduced dynamics (in epsilon)
data.mix = 'R1';
for m = 2:k+1 % zeroth order in R1, no order k red. dyn.
    if ~isempty(R1.R(k-m+2).coeffs)
        W0R1 = W0R1 + coeffs_mixed_terms(k,m, W0,R1.R,data,'R1');
    end
end

WR = W0R1+W1R0;
end
