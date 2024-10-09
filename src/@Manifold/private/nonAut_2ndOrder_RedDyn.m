function [Rk] = nonAut_2ndOrder_RedDyn(I,F,THETA,PHI,Lambda,Lambda_K, M,C,Vm,Ym,l,z_k,order)
% NONAUT_2NDORDER_REDDYN This function computes the reduced dynamics coefficients  
%
% at order k for harmonic i using the second order system routine.
%
% [Rk] = NONAUT_2NDORDER_REDDYN(I,F,THETA,PHI,Lambda,Lambda_K, M,C,Vm,Ym,l,z_k,order)
%
% I:        Eigenvalue positions that are resonant
% F:        Multi-index positions that are resonant for order > 0,
%           otherwise the frequency harmonic
% Lambda:   vector containing master spectrum
% Lambda_K: sum(Lambda * K + kappa * Omega)
% M:        Mass matrix
% C:        Damping Matrix
% Vm:       Velocity part of the lower order and forcing RHS
% Ym:       Displacement part of the lower order and forcing RHS
% l:        SSM dimension
% z_k:      number of multi-indices at order k
% order:    current order of computation
%
% R1:       non-autonomous RD coefficients
%
% See also: NONAUT_2NDORDER_LEADTERMS, NONAUT_2NDORDER_HIGHTERMS

Rk = zeros(l,z_k);

if order == 0
    
    THETA_I = THETA(:,I);
    
    RHS   = -THETA_I' * (Lambda_K*M*Vm +Ym);
    C_0_r = eye(size(RHS,1)); % Fix leading order coefficient to exact resonance    
    Rk(I,:) = lsqminnorm(C_0_r,RHS);
    
else
    % unique multi indices
    if any(F)
        [F_un, ~,i_F_un] = unique(F.');
        
        % Loop over multi-indices that lead to resonance
        fi = 1;
        
        for f = F_un
            I_f = I((i_F_un == fi)); % All resonant eigenvalues for this multi - index
            THETA_f = THETA(:,I_f);
            PHI_f   = PHI(:,I_f);
            
            RHS   = - Lambda_K(f).*THETA_f' *M*Vm(:,f) -THETA_f'*Ym(:,f);
            C_k_r =  THETA_f' * ( C* PHI_f + M * ((Lambda_K(f) + Lambda(I_f).') .* PHI_f)); %Coefficient Matrix
            Rk(I_f,f) = lsqminnorm(C_k_r,RHS);            
            
            % tbd: check implementation for 1:1 internal resonances!           
         
            fi = fi +1;
        end
    end
    
end
end