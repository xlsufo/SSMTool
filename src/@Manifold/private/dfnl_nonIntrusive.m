function [dfnl] = dfnl_nonIntrusive(fun, nlorder,W,X,M,data,nKappas)
% DFNL_NONINTRUSIVE This function computes the nonlinear contributions of internal forces
%
% to the non-autonomous invariance equation for multi indices in set K, given a non-intrusive function handle. 
% This function assumes the displacement parametrisation to be in the first
% n entries of the full SSM parametrisation, ie a first order system of the
% form z = [x ; dot(x)]
%
% [dfnl]= DFNL_NONINTRUSIVE(fun, nlorder,W,X,M,data,nKappas)
%
% fun:      non-intrusive function handle for the jacobian of the nonlinearity
% W:        autonomous SSM coefficients
% X:        non-autonomous SSM coefficients
% nl_order: order of nonlinearity considered
% M:        set of multi indices at current order of computation
% data:     struct containing necessary information for the computation
% nKappas:  number of harmonics present
%
% dfnl:     Jacobian composed with aut SSM parametrisation W acted on non-autonomous coefficients X, contributions at
%           order K
%
% See also: NONAUT_2NDORDER_HIGHTERMS, NONAUT_1STORDER_HIGHTERMS, DFNL_SEMIINTRUSIVE


N = size(W(1).coeffs,1);
dfnl = repmat(struct('val' ,sparse(N,size(M,2))),nKappas, 1);

ordering = data.ordering;

% dimensionality of input vector for nonlinearity
mx_idx = data.nl_input_dim;

% If SSM dimension is uneven, overrule and use full system dimension
if rem(N,2) == 1 
    mx_idx = N;
end

for i =1:size(M,2)
    % All combinations of 2 multi-indices K and L to sum up to any M(:,i)
    
    [KpL,~,~,~] = multi_nsumk(2,M(:,i)); 
    KpL = KpL{1};

    
    DFm = repmat(struct('val',sparse(N,1)),nKappas,1);
    Xjj = repmat(struct('val',sparse(mx_idx,1)),nKappas,1);
    for j = 1:size(KpL,3)
        KpL_abs = sum(KpL(:,:,j));
        
        if KpL_abs(1) == 0
            % Nonlinearity is trivial for trivial input
            continue
        elseif KpL_abs(1) < nlorder-1
            % DF(n) gives minimal nl order of nlorder-1 and SSM coeffs are
            % at least linear.
            continue
        end
        
        K = KpL(:,1,j);
        L = KpL(:,2,j);
        
        L_idx = multi_index_2_ordering(L,ordering,[]);
             
        
        for jj = 1:nKappas
            Xjj(jj).val = X(jj).W(sum(L)+1).coeffs(1:mx_idx,L_idx);
        end
        
        [DFk] = compute_DF(fun,nlorder,W,K,mx_idx,N,ordering,Xjj,nKappas);
        
        for jj = 1:nKappas
            DFm(jj).val = DFm(jj).val + DFk(jj).val;
               
        end
        
        
    end
    
    for jj = 1:nKappas
    dfnl(jj).val(:,i) = DFm(jj).val;
    end
end

end


function [DFk] = compute_DF(dfun,nlorder,W,K,mx_idx,N,ordering,X,nKappas)
% Computes DF ( W ) for every multi-index K(:,m)

% Redundantly also computes combos with zero vectors
[g,~,~,~] = multi_nsumk(nlorder-1,K,'unique');
g = g{1};

DFk = repmat(struct('val',sparse(N,1)),nKappas,1);

% Loop over all partitions of m
for i = 1:size(g,3)    
    h_abs = sum(g(:,:,i));
    gi = g(:,:,i);
    % Take out combinations which have a zero vector
    if any(h_abs == 0)
        continue
    else
        h_idx = multi_index_2_ordering(g(:,:,i),ordering,[]);
        
        vectors = cell(nlorder-1,1);
        
        % Check for all zero SSM coefficient vector
        emptyflag = false; 
        for j = 1:(nlorder-1)

            if ~isempty(W(h_abs(j)).coeffs)
                vectors{j} = W(h_abs(j)).coeffs(1:mx_idx,h_idx(j));  
            else
                emptyflag = true;
                continue
            end
            
        end
        
        if emptyflag %no contribution of nl
            continue
        end

        for jj = 1:nKappas
            
            DF_comp_W =  StEP(dfun, vectors, gi, nlorder-1 );            
            DFk(jj).val = DFk(jj).val + DF_comp_W(:,1:mx_idx)*X(jj).val;
            
        end
    end
end
end

