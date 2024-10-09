function [sigma_out] = dfnl_intrusive(N_ind,W1 ,m,data,nKappas)
% DFNL_INTRUSIVE This function computes the nonlinear contributions to the 
%
% non-autonomous invariance equation for multi indices in set K, given a intrusive tensor fnl. 
% The Sigmas are contributions to the jacobian of the nonlinearity. Details
% can be found in  $\textit{Computing SSMs using multi indices - nonautonomous}$
%
% [sigma_out] = DFNL_INTRUSIVE(N_ind,W1 ,m,data,nKappas)
%
% N_ind:    multi indices of function which is being composed with
%           parametrisation, the jacobian of fnl
% W1:       non-autonomous SSM coefficients
% m:        order of computation
% data:     data struct containing necessary information for computation
% nKappas:  number of harmonics in forcing
%
% sigma_out:composed jacobian with autonomous SSM and acted on
%           non-autonomous SSM, at order m, for all harmonics
%           
% See also: NONAUT_FEXT_PLUS_FNL, FNL_INTRUSIVE, MULTI_NSUMK


%% Vectorised computation

M = flip(sortrows(nsumk(data.l,m,'nonnegative')).',2);

sigma_out = repmat(struct('val',zeros(size(N_ind,2),size(M,2))),nKappas,1);

N_abs = sum(N_ind);
One = speye(data.N);

<<<<<<< HEAD:src/@Manifold/private/dfnl_intrusive.m
[NOne_base,idx_N_base,idx_One_Base] = multi_subtraction(N_ind,One,'Identity');
for k = 0:m-1 % starting at 0, order k leads to zero pi
    K = W1(1).W(k+1).ind;
    if ~isempty(K)
        cond = (m-(sum(K(:,1)))) >= N_abs -1;%otherwise pi is zero anyways
=======
N_abs = sum(N);
One = speye(field.N);

[NOne_base,idx_N_base,idx_One_Base] = multi_subtraction(N,One,'Physical');
for m = 0:k-1 % starting at 0, order k leads to zero pi
    M = W(m+1).ind;
    if ~isempty(M)
        cond = (k-(sum(M(:,1)))) >= N_abs -1;%otherwise pi is zero anyways
>>>>>>> main:src/@Manifold/private/compute_sigma.m
        if nnz(cond)>0
            [L,idx_K,idx_M] = multi_subtraction(M,K,'Arbitrary');
            
            %find which multi indices in NOne correspond to multi indices
            %in N that fulfill cond
            idx_base = ismember(idx_N_base,find(cond));
            [NOne,idx_N,idx_One] = deal(NOne_base(:,idx_base),idx_N_base(idx_base),idx_One_Base(idx_base));
            
            pis = fnl_intrusive(NOne,L,data);
            
            % Sort columns in ascending order of index in K
            [idx_K,sort_idx_K] = sort(idx_K);
            idx_M = idx_M(sort_idx_K);
            pis = pis(:,sort_idx_K);
            
            % Sort rows in ascending order of index in N
            [idx_N,sort_idx_N] = sort(idx_N);
            idx_One = idx_One(sort_idx_N);
            pis = pis(sort_idx_N,:);
            
            % Read out corresponding SSM coefficients and multiply by pi
            % (redundancy could be resolved by sorting by idx_M and only reading out each once )
            
            [X_w,Y_w] = meshgrid(idx_One,idx_M);
            X_w = reshape(X_w,[],1);
            Y_w = reshape(Y_w,[],1);
            lin_idx_W = sub2ind(size(W1(1).W(k+1).coeffs),X_w,Y_w);
            
            %% Results for each harmonic
            for jj = 1:nKappas
                            
            W_coeffs = W1(jj).W(k+1).coeffs(lin_idx_W);            
            W_pi = full(reshape(W_coeffs,numel(idx_K),numel(idx_N)).'.* pis);

            
            % Sum all the contributions that correspond to the same columns
            % in K
            [xx, yy] = ndgrid(idx_K.',1:size(W_pi, 1));
            accum_Wpi = accumarray([yy(:) xx(:) ], reshape(W_pi.', 1, []));
                           
            % Getting all (n_j)s
            % (could be more efficient if redundancy in n_js is removed)
            lin_idx_N = sub2ind(size(N_ind),idx_One,idx_N);
            accum_Wpi = accum_Wpi .*reshape(N_ind(lin_idx_N),[],1);
                        
            
            %nposition of each element in accum_Wpi wrt N(:,cond) for
            %summing them
            [~,ia,ic] = unique(idx_N,'stable');
            pos_idx_N = 1:numel(ia);
            pos_idx_N = pos_idx_N(ic); 
            
            if size(pos_idx_N,2)>1
                pos_idx_N = pos_idx_N.';
            end
            
            % Summing terms for same column in N
            [xx, yy] = ndgrid(pos_idx_N,1:size(accum_Wpi, 2));
            sigma_m = accumarray([xx(:) yy(:)], accum_Wpi(:));
                       
            %Create linear indices
            [X,Y] = meshgrid(unique(idx_N),unique(idx_K));
            X = reshape(X,[],1);           
            Y = reshape(Y,[],1);
            lin_idx = sub2ind(size(sigma_out(1).val),X,Y);

            % Read out results
            
            if size(N_ind,2) > 1
                sigma_out(jj).val(lin_idx) = sigma_out(jj).val(lin_idx) + reshape(sigma_m.',[],1);
            else
                sigma_out(jj).val(lin_idx) = sigma_out(jj).val(lin_idx) + reshape(sigma_m.',1,[]);
            end
            end
        end
    end

end

end

