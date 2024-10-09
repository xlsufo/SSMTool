function [R1_ik,W1_ik] = nonAut_1stOrder_SolveInvEq(FG,WR,data)
% NONAUT_1STORDER_SOLVEINVEQ This function sets up the non-autonomous invariance  
%
% equation for the first order system computation and solves it with the
% predefined solver.
%
% [R1_ik,W1_ik] = NONAUT_1STORDER_SOLVEINVEQ(FG,WR,data)
%
% FG:       internal and external forces composed with SSM coefficients for
%           harmonic i.
% WR:       Composed terms W1R0 + W0R1 at order k and for harmonic i
% data:     data struct containing necessary information for computation
%
% W1_ik:    non-autonomous SSM coefficients at order k for harmonic i
% R1_ik:    non-autonomous RD coefficients at order k for harmonic i
%
% See also: NONAUT_2NDORDER_SOLVEINVEQ, NONAUT_1STORDER_HIGHTERMS

% Unpack variables
kappas = data.kappas;
Omega  = data.Omega;
W_M    = data.W_M;
A      = data.A;
B      = data.B;
i      = data.i;
k      = data.k;
I      = data.I;
F      = data.F;
l      = data.l;
N      = data.N;
K_lambda = data.Lambda_K;
W0     = data.W0;
solver = data.solver;

z_k = nchoosek(k+l-1,l-1);

% Computes the SSM coeffs and reduced dynamics coefficients using first
% order implementation

R1_ik_coeff     = sum( conj(W_M(:,I)).* ( FG(:,F) - B*(WR(:,F))));
R1_ik           = sparse(I,F,R1_ik_coeff , l,nchoosek(k+l-1,l-1));
R1_dum(k+1).coeffs    = R1_ik;

if l > 1
    R1_dum(k+1).ind    = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
else
    R1_dum(k+1).ind = k;
end

%% Solve the coefficient equation for the SSM coefficients
% Add R1 order k contribution to the right hand side

RHS                 = B* (WR +  coeffs_mixed_terms(k,1, W0,R1_dum,data,'R1')) - FG;

W1_ik = zeros(N,z_k);

for j = 1:z_k
    if any(RHS(:,j))
        C_i = A - B * (K_lambda(j) + 1i * kappas(:,i)*Omega); % Coefficient matrix
        W1_ik(:,j) = solveinveq(C_i,RHS(:,j),solver);
    end
end
end