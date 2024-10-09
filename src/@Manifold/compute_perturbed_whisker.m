function [W1, R1, varargout] = compute_perturbed_whisker(obj, order,W0,R0,varargin)
% COMPUTE_PERTURBED_WHISKER This function computes the non-autonomous SSM   
%
% up to order 'order'.
%
% [W1, R1] = COMPUTE_PERTURBED_WHISKER(obj, order,W0,R0,varargin)
%
% obj:      SSM class object
% order:    approximation order up until which SSM is computed
% W0:       autonomous SSM coefficients
% R0:       autonomous RD coefficients
% varargin: if not stored in the obj class, this can be used to input Omega
%
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
%
% See also: COMPUTE_WHISKER, NONAUT_1STORDER_WHISKER, NONAUT_2NDORDER_WHISKER

%% Non-autonomous (quasi)periodic perturbation to whiskers of invariant manifolds
% We consider the mechanical system
%
% $$$\mathbf{B}\dot{\mathbf{x}} =\mathbf{Ax}+\mathbf{G}(\mathbf{x})+\epsilon\mathbf{F}(\mathbf{\phi},
% \mathbf{x}),$$$$\dot{\mathbf{\phi}}	=\mathbf{\Omega}$$
%
% with quasi-periodic forcing.
%
% In the non-autonomous setting, the SSM and the corresponding reduced dynamics
% would be parameterized by the angular variables $\mathbf{\phi}$, as well. In
% general, we may write
%
% $$    \mathbf{S(p,{\mathbf{\phi}})} = \mathbf{T}\mathbf{(p)} + \epsilon \mathbf{U}\mathbf{(p,\mathbf{\phi})}
% + O(\epsilon^2),$$
%
% $$    \mathbf{R(p,{\mathbf{\phi}})} = \mathbf{P}\mathbf{(p)} + \epsilon \mathbf{Q}\mathbf{(p,\mathbf{\phi})}
% + O(\epsilon^2),$$
%
% where $\mathbf{T}(\mathbf{p}),\mathbf{P}(\mathbf{p})$ recover the SSM and
% reduced dynamics coefficients in the unforced limit of $\epsilon=0.$
%
% These functions as well as the nonlinearity and the forcing are expanded in
% phase space coordinates. The time dependent coefficients of those expansions
% are furthermore expanded as a Fourier-series. As an example, the Force and the
% non-autonomous SSM-coefficients are given as
%
% $$    \mathbf {F}(\mathbf{x},\mathbf{\phi}) =     \left[  f^1(\mathbf{x},\mathbf{\phi}),
% \cdots ,   f^{2n}(\mathbf{x},\mathbf{\phi})     \right]^T,     \ f^i(\mathbf{x},\mathbf{\phi})
% = \sum_{\mathbf{n}\in \mathbb{N}^{2n}} F^i_{\mathbf{n}}(\mathbf{\phi}) \mathbf{x}^\mathbf{n}$$
%
% $$F^b_{\mathbf{k}}(\mathbf{\phi}) = \sum_{\mathbf{\eta} \in \mathbb{Z}^k}
% F^b_{\mathbf{k},\mathbf{\eta} } e^{i\langle \mathbf{\eta},\mathbf{\phi}\rangle}$$
%
% $$    \mathbf {U}(\mathbf{p},\mathbf{\phi}) =     \left[  u^1(\mathbf{p},\mathbf{\phi}),
% \cdots ,   u^{2n}(\mathbf{p},\mathbf{\phi})     \right]^T,     \ u^i(\mathbf{p},\mathbf{\phi})
% = \sum_{\mathbf{m}\in \mathbb{N}^{l}} U^i_{\mathbf{m}}(\mathbf{\phi}) \mathbf{p}^\mathbf{m}$$
%
% $$U^i_{\mathbf{k}}(\mathbf{\phi}) = \sum_{\mathbf{\eta} \in \mathbb{Z}^k}
% U^i_{\mathbf{k},\mathbf{\eta} } e^{i\langle \mathbf{\eta},\mathbf{\phi}\rangle}$$
%
% This leads to the invariane equation
%
% $$ \mathbf{B}  \bigg( \text{D}_\mathbf{p}( \mathbf{T}\mathbf{(p)})\mathbf{Q}\mathbf{(p,\mathbf{\phi})}
% + (\partial_\mathbf{p} \mathbf{U}\mathbf{(p,\mathbf{\phi)})  \mathbf{P}\mathbf{(p)}+
% (\partial_\mathbf{\phi} \mathbf{S(p,\mathbf{\phi})})   \mathbf{\Omega }  \bigg)=\\
% \ \ \ \ \mathbf{A}\mathbf{U}\mathbf{(p,\mathbf{\phi})} +     \big[\text{D}_\mathbf{x}\mathbf{G}
% \circ \mathbf{T}(\mathbf{p}) \big]\mathbf{U}(\mathbf{p},\mathbf{\phi})+
% \mathbf{F} (\mathbf{\phi},\mathbf{S(p,\mathbf{\phi})})$$
%
% The various expansions are plugged into this equation and then the equation
% is iteratively solved for the coefficients. The functions $\mathbf{T(p),P(p)}$
% are already known from the autonomous computation, their coefficients given
% by $\texttt{W0}$ (SSM) and $\texttt{R0}$ (reduced dynamics) .
%
% The external force $\mathbf{F}$ is input as a field of the property $\texttt{System}$
% of the SSM object. Since the equations for different frequency multi-indices
% decouple and the code is parallelised over these decoupled equations we want
% to make the read out hirarchy such that the first parameter corresponds to the
% frequency multi-indices. $\texttt{System.Fext.data(i)}$ indices into the $\texttt{i}$-th
% component of the field $\texttt{data}$ which is a struct array containing struct
% arrays. There is one such  contained struct array for each frequency multi-index.
%
% Every struct now contains two arrays with the coefficients and the spatial
% multi-indices respectively. The coefficients of the force for the $i$-th frequency-
% and ther order $k$ spatial multi-indices and their coefficients are stored in
% the rows of  $\texttt{data(i).F\_n\_k(k).ind} $ and the columns of $\texttt{data(i).F\_n\_k(k).coeffs}$
% . The $i$-th frequency multi-index is stored in $\texttt{data(i).kappa}$.
%
% The non-autonomous SSM and reduced dynamics coefficients are stored analogously.
% In $\texttt{W\_1(i).W(k).coeffs}$ the coefficients of the SSM expansion corresponding
% to $\mathbf{\eta}_i$ and order $k$ spatial multi-indices are stored. During
% the computation the multi-indices are stored in the columns of $\texttt{W\_1(i).W(k).ind}$
% in reverse lexicographic order, upon outputting the resulting coefficients however
% the storing scheme is reversed, in the output the multi-indices are stored in
% the rows in lexicographic ordering, the standard way of storing used throughout
% the software package.
%
% While in the documentation the frequency multi-indices are called $\mathbf{\eta}$
% for good distinguishability from spatial multi-indices in the code they are
% called $\texttt{kappa}$.
%% System Properties


if isempty(varargin)
    Omega = obj.System.Omega; 
else
    Omega = varargin{1};
end

data.Omega  = Omega;
data.order  = order;
data.nKappa = obj.System.nKappa;
data.l      = obj.dimManifold;        % dimension of manifold
data.N      = obj.dimSystem;          % phase space size
data.ordering = 'revlex';
data.Lambda_M_vector = obj.E.spectrum;
data.solver = obj.Options.solver;
data.reltol = obj.Options.reltol;

if ~obj.Options.contribNonAuto 
    obj.solInfoNonAut.timeEstimate = 0;
    
elseif isempty(obj.solInfoNonAut.timeEstimate) && order > 0
    
    % struct to save computation times
    obj.solInfoNonAut.timeEstimate = zeros(1,order+1);
    obj.solInfoNonAut.nlTime  = zeros(1,order+1);
    obj.solInfoNonAut.mixTime = zeros(1,order+1);
    obj.solInfoNonAut.eqTime  = zeros(1,order+1);    
end

if strcmp(obj.Options.COMPtype,'first')
    % First order Computation
    [W1,R1,nRHS] = nonAut_1stOrder_whisker(obj,W0,R0,data);
    varargout{1} = nRHS;
else
    % Second order computation
    [W1,R1] = nonAut_2ndOrder_whisker(obj,W0,R0,data); 
end

% As these coefficients depend explicitly on omega at higher orders save it
% as property
W1(1).Omega = Omega;
R1(1).Omega = Omega;

% Set correct ordering of coefficients
if obj.Options.contribNonAuto
    compOrder = order;
else
    compOrder = 0;
end

for i = 1:data.nKappa    
    % Output coefficients in lexicographic ordering, with multi indices stored
    % in rows
    for k = 1:compOrder+1 %index starts at 0
        W1(i).W(k) = coeffs_lex2revlex(W1(i).W(k),'TaylorCoeff');
        R1(i).R(k) = coeffs_lex2revlex(R1(i).R(k),'TaylorCoeff');
    end
end

end
