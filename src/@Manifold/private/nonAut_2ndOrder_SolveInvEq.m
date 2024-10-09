function [R1_ik,W1_ik] = nonAut_2ndOrder_SolveInvEq(FG,WR,data)
% NONAUT_2NDORDER_SOLVEINVEQ This function sets up the non-autonomous invariance  
%
% equation for the second order system computation and solves it with the
% predefined solver.
%
% [R1_ik,W1_ik] = NONAUT_2NDORDER_SOLVEINVEQ(FG,WR,data)
%
% FG:       internal and external forces composed with SSM coefficients for
%           harmonic i.
% WR:       Composed terms W1R0 + W0R1 at order k and for harmonic i
% data:     data struct containing necessary information for computation
%
% W1_ik:    non-autonomous SSM coefficients at order k for harmonic i
% R1_ik:    non-autonomous RD coefficients at order k for harmonic i
%
% See also: NONAUT_1STORDER_SOLVEINVEQ, NONAUT_2NDORDER_HIGHTERMS

% Unpack data
kappas = data.kappas;
Omega  = data.Omega;
[i,k,I,F,THETA,PHI,Lambda_K,Lambda_M_vector,Mass,Damp,Stiff,solver,l,N] = deal(data.i,data.k,data.I,data.F,data.THETA,data.PHI, ...
    data.Lambda_K,data.Lambda_M_vector,data.Mass,data.Damp,data.Stiff,data.solver,data.l,data.N);
n = N/2;

%% Start computation
Ym     = ( Damp * WR(1:n,:) + Mass * WR((n+1):end,:)) - FG(1:n,:);
Vm     = WR(1:n,:);

Lambda_K_Om = Lambda_K + 1i * (kappas(:,i)*Omega);

z_k = nchoosek(k+l-1,l-1);

[Rk] = nonAut_2ndOrder_RedDyn(I,F,THETA,PHI,Lambda_M_vector,Lambda_K_Om, Mass,Damp,Vm,Ym,l,z_k,k);

w_1i    = zeros(n,z_k);
w_1idot = zeros(n,z_k);

for f = 1:z_k
    
    L_k = ( Mass * ((Lambda_K_Om(f) + Lambda_M_vector.') .* PHI) + Damp*PHI ) * Rk(:,f);
    L_k = L_k + Lambda_K_Om(f)*Mass*Vm(:,f) + Ym(:,f);
    
    C_k = -(Stiff + Lambda_K_Om(f)*Damp + Lambda_K_Om(f)^2 *Mass );
    
    if any(L_k)
        w_1i(:,f)    = solveinveq(C_k,L_k,solver);        
    end
    
    w_1idot(:,f) = Lambda_K_Om(f) * w_1i(:,f) + PHI * Rk(:,f) + Vm(:,f);
end

W1_ik       = [w_1i;w_1idot];
R1_ik       = sparse(Rk);

end
