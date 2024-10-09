function [E, I_k,K_lambda]       = nonAut_resonant_terms(k,kappa,data,order)
% NONAUT_RESONANT_TERMS  This function finds the combinations of frequency multi-indices,
%
% master mode eigenvalues and the spatial multi-indices at zeroth and order 
% k that lead to internal resonances.
%
% [E, I_k,K_lambda]       = NONAUT_RESONANT_TERMS(k,kappa,data,order)
%
% k:        current order of SSM computation
% kappa:    at zeroth order - vector containing all harmonics
%           at higher order - the harmonic for which SSM coefficients are currently computed
% data:     data struct containing necessary information for computation
% order:    approximation order up until which SSM is computed
%
% E:        indices of the master modes that lead to resonances        
% I_k:      at zeroth order  - index of kappa that leads to resonance
%           at higher orders - index of multi-indices that lead to resonance
% K_lambda: direct product of multi-indices at order k with master mode
%           vector
%
% See also: AUT_RESONANT_TERMS 

Omega  = data.Omega;
Lambda = data.Lambda_M_vector;   % master modes eigenvalues
l      = data.l;
% Tolerance for resonances
ref = min(abs(Lambda));
abstol = data.reltol * ref;

switch order
    case 'zero'
        %% Find zeroth order resonant terms
        % We determine the near inner resonances of the coefficient matrix where
        %
        % $$     \lambda_{j} - i\langle \mathbf{\eta}_{f}, \mathbf{\Omega } \rangle\approx
        % 0, \ e\in \{1,...,l\}, \ f \in\{1,...,K\}$$
        %
        % holds. The index pairs that fulfill this condition are stored.
        %
        % $$E := \{ e_1,  ... ,e_{r_{ext}} \in \{1,...,l\}\} \\ F := \{ f_1,  ... ,f_{r_{ext}}
        % \in \{1,...,K\}\}$$
        
        % kappa in this case contains all kappas
        lambda_C_10 =  repmat(Lambda,[1,size(kappa,2)]) - 1i*repmat(kappa*Omega,[l 1]);
        
        [E, I_k] = find(abs(lambda_C_10)<abstol);
        K_lambda = [];
        
        % I_k contains the frequency index
    case 'k'
        %% Find higher order resonant terms
        % The coefficient matrix for frequency multi-index $\mathbf{\eta}$ shows singularities
        % if the resonance condition
        %
        % $$    \lambda_e - \bigg( \sum_{j=1}^l k_j\lambda_j     + i \langle \mathbf{\Omega},
        % \mathbf{\eta} \rangle \bigg) \approx 0$$
        %
        % is fulfilled for some $\lambda_e$ in the master subspace. We therefore have
        % to find all such resonant combinations.
        
        %Find the resonances
        
        if l > 1
            K = flip(sortrows(nsumk(l,k,'nonnegative')).',2); %order k multi-indices
        else
            K = k;
        end
        z_k = size(K,2);
        %vector with each element korresponding to summing multi_index k with all master lambdas
        K_lambda = sum(K .* Lambda);
        lambda_C_11 = repmat(Lambda,[1,z_k]) - repmat(K_lambda + 1i * (kappa*Omega),[l 1]);
        
        [E, I_k] = find(abs(lambda_C_11)<abstol); %I_k indicates the spatial multi-index the resonance corresponds to
end
end