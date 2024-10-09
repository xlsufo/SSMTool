function [W1,R1,varargout] = nonAut_1stOrder_leadTerms(W1,R1,data)
% NONAUT_1STORDER_LEADTERMS This function computes the leading order,
%
% parametrisation-coordinate independent parts of the non-autonomous SSM 
% using the second order system computation routine.
%
% [W1,R1, varagout] = NONAUT_1STORDER_LEADTERMS(W1,R1, data)
%
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
% data:     data struct containing necessary information for computation
%
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
% varargout:information concering the solution time of the inv. equation.
%
% See also: NONAUT_2NDORDER_LEADTERMS, NONAUT_1STORDER_WHISKER

% Unpack variables
Omega  = data.Omega;
nKappa = data.nKappa;
N      = data.N; 
l      = data.l;   
kappas = data.kappas;

% System matrices
A     = data.A;      
B     = data.B;    
FEXT  = data.FEXT;
% Manifold variables

W_M   = data.W_M;
V_M   = data.V_M;


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

%% Set reduced dynamics
% These sets now determine the bases of the near left kernel of the coefficient
% matrices $\mathbf{\mathcal{L}_{0,\eta}\cdot$ onto which we project the RHS of
% the coefficient equation and set it equal zero. This gives the explicit expression
% for the reduced dynamics coefficients.
%
% $$    Q^{e_i}_{\mathbf{0},\mathbf{\eta}_{f_i}}    =    \langle     \mathbf{w}_{e_i}
% ,     	\big[     		F^1_{\mathbf{0},\mathbf{\eta}_{f_i}}      		     		\cdots
% F^{2n}_{\mathbf{0},\mathbf{\eta}_{f_i}}     	\big]^T      \rangle$$

if r_ext
    Q_0 = sparse(ev_idx,harm_idx,sum(conj(W_M(:,ev_idx)).* F_0(:,idx_0(harm_idx))), l,numel(idx_0));
    RHS = B*V_M*Q_0 - F_0(:,idx_0);
else
    Q_0 = sparse(l,nKappa);
    RHS = - F_0(:,idx_0);
end
varargout{2} = -RHS; % to be used in optimization
%% Solve for the order zeroth order SSM-coefficients W1, R1

run_idx = 1;
for j= idx_0
    
    R1(j).R(1).coeffs = Q_0(:,run_idx);
    R1(j).R(1).ind    = sparse(l,1);
    run_idx = run_idx + 1;
end

soltic = tic;
if data.NonAuto % whether to ignore higher order
    kappas_0 = kappas(:,idx_0);
    [redConj,mapConj] = nonAut_conj_red(kappas_0, F_0(:,idx_0));
    
    for j = 1:numel(redConj)
        
        C_j  =  A - 1i*dot(Omega,kappas_0(:,redConj(j)))*B ;
        
        W10j = solveinveq(C_j,RHS(:,redConj(j)),data.solver);
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
varargout{1}= soltime;
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

