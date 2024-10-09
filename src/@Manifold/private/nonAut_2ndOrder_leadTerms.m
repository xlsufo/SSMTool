function [W1,R1,varargout] = nonAut_2ndOrder_leadTerms(W1,R1,data)
% NONAUT_2NDORDER_LEADTERMS This function computes the leading order,
%
% parametrisation-coordinate independent parts of the non-autonomous SSM 
% using the second order system computation routine.
%
% [W1,R1, varagout] = NONAUT_2NDORDER_LEADTERMS(W1,R1, data)
%
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
% data:     data struct containing necessary information for computation
%
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
% varargout:information concering the solution time of the inv. equation.
%
% See also: NONAUT_1STORDER_LEADTERMS, NONAUT_2NDORDER_WHISKER

% Unpacking variables

Omega  = data.Omega;
nKappa = data.nKappa;
N      = data.N; 
n      = N/2;
l      = data.l;   
kappas = data.kappas;

% System matrices
Mass  = data.Mass;      
Damp  = data.Damp;    
Stiff = data.Stiff;
FEXT  = data.FEXT;
% Manifold variables
THETA           = data.THETA;
PHI             = data.PHI;
Lambda_M_vector = data.Lambda_M_vector;

%% Solving for coefficients with k=0
% The coefficient equation for this case reads
%
% $$\sum_{i=1}^{2n}     \underbrace{    \big(        \mathbf{(A)}_{bi}  - i
% \langle \mathbf{\eta},\mathbf{\Omega }\rangle        \mathbf{B}_{bi}    \big)
% }_{:= (\mathcal{L}_{\mathbf{0},\mathbf{\eta})_{bi}}}    U^i_{\mathbf{0},\mathbf{\eta}}=\sum_{j=1}^{l}
% (\mathbf{Bv}_j)_b      Q^j_{\mathbf{0},\mathbf{\eta}}-    F_{b,\mathbf{0},\mathbf{\eta}}$$

% Finding all force contributions at zeroth order
[F_0, idx_0] = zeroth_order_forcing(nKappa,N,FEXT);
%% Find resonant terms

[ev_idx, harm_idx,~] = nonAut_resonant_terms([],kappas(:,idx_0),data,'zero'); % contains harmonic idx.
r_ext    = length(ev_idx);

V0     = sparse(n,1);

if r_ext
    % unique harmonic positions
    [harm_idx_un, ~,i_E_un] = unique(harm_idx.'); 
    
    if numel(harm_idx_un) > 0
        Y0 = zeros(N/2,numel(idx_0));
    else
        Y0 = sparse(N/2,numel(idx_0));
    end
    
    jj = 1;
    for un_harm = harm_idx_un % Loop over all unique resonant harmonics
        un_ev = ev_idx((i_E_un == jj)); % All resonant eigenvalues for this harmonic
        
        Y0_tmp     = - F_0(1:n,un_harm);
        Y0(:,un_harm) = Y0_tmp;
        Lambda_0_Om =  1i * (kappas(:,idx_0(un_harm))*Omega);
        
        [S0] = nonAut_2ndOrder_RedDyn(un_ev,[],THETA,PHI,Lambda_M_vector,Lambda_0_Om, Mass,Damp,V0,Y0_tmp,l,1,0);
        
        R1(idx_0(un_harm)).R(1).coeffs = sparse(S0);
        R1(idx_0(un_harm)).R(1).ind    = sparse(l,1);
        jj = jj +1;
    end        
end

soltic = tic;
if data.NonAuto % whether to ignore higher order
    kappas_0 = kappas(:,idx_0);
    [redConj,mapConj] = nonAut_conj_red(kappas_0, F_0(:,idx_0));
    
    for j = 1:numel(redConj)
        
        
        S0 = R1(j).R(1).coeffs;
        if isempty(S0)
            S0 = 0;
        end
        
        Lambda_0_Om =  1i * (kappas(:,redConj(j))*Omega);
        
        L_k = ( Mass * ((Lambda_0_Om + Lambda_M_vector.') .* PHI) + Damp*PHI ) * S0;
        L_k = L_k +  Y0(:,redConj(j));
        
        C_k = -(Stiff + Lambda_0_Om*Damp + Lambda_0_Om^2 *Mass );
        w_10    = solveinveq(C_k,L_k,data.solver);
        w_10dot = Lambda_0_Om * w_10 + PHI * S0 + V0;
        
        W10j = [w_10;w_10dot];
        mapj = mapConj{j};
        
        
        switch numel(mapj)
            case 1
                W1(idx_0(mapj(1))).W(1).coeffs = W10j;
                W1(idx_0(mapj(1))).W(1).ind    = sparse(l,1);
            case 2
                
                W1(idx_0(mapj(1))).W(1).coeffs = W10j;
                W1(idx_0(mapj(2))).W(1).coeffs = conj(W10j);
                
                W1(idx_0(mapj(1))).W(1).ind    = sparse(l,1);
                W1(idx_0(mapj(2))).W(1).ind    = sparse(l,1);
            otherwise
                error('there exist redundancy in kappa of external forcing');
        end
    end
    
end
soltime = toc(soltic);
varargout{1} = soltime;
end

function [F_0, idx_0] = zeroth_order_forcing(nKappa,N,FEXT)
% Finding all force contributions at zeroth order
F_0 = zeros(N,nKappa);
for i = 1:nKappa
    if ~isempty(FEXT.data(i).F_n_k(1).coeffs)
        F_0(:,i)= FEXT.data(i).F_n_k(1).coeffs;  % each column corresponds to one kappa
    end
end
idx_0 = find(any(F_0~=0)); % index for all kappas that contribute at zeroth order
end
