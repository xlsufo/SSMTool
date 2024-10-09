function [Rk] = Aut_2ndOrder_RedDyn(I,F,THETA,PHI,C,Lambda_K,Lambda, M,Vm,Ym,l,z_k)
% AUT_2NDORDER_REDDYN
%
% This function computes the autonomous reduced dynamics coefficients at
% order k using the second order systems routine. 
%
% [Rk] = AUT_2NDORDER_REDDYN(I,F,THETA,PHI,C,Lambda_K,Lambda, M,Vm,Ym,l,z_k)
%
% F:        multi-index positions that are resonant
% I:        Eigenvalue positions that are resonant
% THETA:    Left Eigenvectors of system
% PHI:      Right Eigenvectors of system
% C:        Damping matrix of system
% Lambda_K:  direct product of multi-indices at order K with master mode
%           vector
% Lambda:   Vector containing master mode eigenvalues
% M:        Mass matrix
% Vm:       Velocity variable RHS of invariance equation
% Ym:       Displacement variable RHS of invariance equation
% l:        dimension of SSM
% z_k:      number of distinct multi-indices at order k
%
% Rk:       order k autonomous reduced dynamics
%
% See also: AUT_2NDORDER_SSM, NONAUT_2NDORDER_REDDYN

Rk = zeros(l,z_k);
% unique multi indices
if any(F)
    [F_un, ~,i_F_un] = unique(F.');
    
    % Loop over multi-indices that lead to resonance
    ii = 1;

    % group the multi-indices by eigenvalues
    % if different eigenvalues resonant for same multi-index: 1:1 resonance

    for f = F_un

        I_f = I((i_F_un == ii)); % All resonant eigenvalues for this multi - index

        if numel(I_f) == 1

            THETA_f = THETA(:,I_f);
            PHI_f   = PHI(:,I_f);
            LHS = THETA_f' * ( C* PHI_f + M * ((Lambda_K(f) + Lambda(I_f).') .* PHI_f));
            RHS =  Lambda_K(f).*THETA_f' * M*Vm(:,f) + THETA_f' *(Ym(:,f));
            Rk(I_f,f) = -(RHS) / (LHS);
        
        elseif numel(I_f) == 2
            %% Case of 1:1 resonance

            error("Systems with 1:1 internal resonance are not yet supported for 'COMPtype=second'")
            % dummy matrices

            theta_1 = THETA(:,I_f(1));
            theta_2 = THETA(:,I_f(2));
            lambda_1 = Lambda(I_f(1));
            lambda_2 = Lambda(I_f(2));
            phi_1   = PHI(:,I_f(1));
            phi_2   = PHI(:,I_f(2));

            A_11 = theta_1' * (C + M * (Lambda_K(f) + lambda_1)) * phi_1;

            A_12 = theta_1' * (C + M * (Lambda_K(f) + lambda_2)) * phi_2;

            A_21 = theta_2' * (C + M * (Lambda_K(f) + lambda_1)) * phi_1;

            A_22 = theta_2' * (C + M * (Lambda_K(f) + lambda_2)) * phi_2;

            LHS = [A_11 A_12 ; A_21 A_22];
            RHS =  Lambda_K(f).*THETA(:,I_f)' * M*Vm(:,f) + THETA(:,I_f)' *(Ym(:,f));

            idx = sub2ind(size(Rk),I_f,[f;f]);

            Rk(idx) = - LHS \ RHS;

        else
            error("1:1:1 internal resonance not yet implemented for COMPtype = second")

        end

        ii = ii +1;

    end
end
end



%% Set analogous to first order case
% Rk(I_f,f)= Lambda(I_f) .* -THETA_f'*M*Vm(:,f) + -THETA_f' *Ym(:,f);
        