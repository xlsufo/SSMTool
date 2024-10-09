function [g, g_comb,nv_un,nv_ic,varargout] = multi_nsumk(nv,v,varargin)
% MULTI_NSUMK  This function finds all combinations of $nv$ positive integer
%
% vectors adding up to another positive integer vector $\mathbf{v}$:
% 
% Let $\mathbf{v} \in \mathbb{N}^l$. It finds all combinations consisting of 
% $n$ positive, nonzero integer vectors such that they sum up to
% 
% $\mathbf{v}$: $\sum_{i=1}^n \mathbf{g}_i = \mathbf{v}$. All of those combinations 
% are stored as 2D arrays, where the individual columns are filled with thevectors, 
% in column $i$ the vector $\mathbf{g}_i$ is stored. The various combinations 
% are stacked upon each other in a 3D array, which contains a third dimension 
% with size equal to amount of combinations of vectors that sum up to $\mathbf{v}$.
% 
% The result is then stored in a cell array, that contains these combinations 
% for all tuples 
% 
% $i,j$ , corresponding to $\texttt{v(:,i), nv(j)}$.
% there can be multiple vectors in v
%
% This functionality is needed to construct the forcing contribution
% coefficients pi_nk in fnl_intrusive. 
%
% [g, g_comb,nv_un,nv_ic,varargout] = MULTI_NSUMK(nv,v,varargin)
%
% nv:       number of integer vectors that should sum up to v
%           can be an array if v is a matrix, then the algorithm performs
%           the partition for each vector in v..
% v:        vector which is partitioned into nv vectors, or matrix wich
%           contains vectors, for each of which the partitions are computed
% varagin:  'rmzero' can be used to remove combinations that contain zero vectors
%           'unique' can be used to remove permutations of combinations
%
% g:        cell array containing all the combos of nv vectors summing up
%           to v. If multiple input vectors are given, each of its z_k rows
%           contains permutations for one vector in v, each column (j) corresponds
%           to nv(j) combinations that sum up to said vector. 
% g_comb:   contains number of combinations that exist for a said tuple i,j
% nv_un:    unique number of combinations that are requested.
% nv_ic:    map nv_un back to positions in nv via this index array.
% varargout:contains number of occurences of a permutation of combinations 
%           for each tuple i,j
%
% See also: MULTI_ADDITION, MULTI_INDEX_2_ORDERING, MULTI_SUBTRACTION, FNL_INTRUSIVE

%%
z_k = size(v,2);
l = size(v,1);
S = cell(1,l);
[nv_un,~,nv_ic] = unique(nv);

g = cell(z_k,size(nv_un,2));
g_comb = zeros(z_k,size(nv_un,2));
n_run = 1;
for n = nv_un
    for f = 1:z_k
        
        for i = 1:l
            %special case
            if v(i,f) == 1 && n == 1
                S{i} = 1;
            else
                S{i} = nsumk(n,v(i,f),'nonnegative').' ;
                %finds all combinations for each index in v
            end
        end
%% 
% $\texttt{combvec}$ finds all possible combinations of the vectors in $\texttt{S}$.
        g_fn = combvec(S{:});
        g_fn = reshape(g_fn,n,l,[]);
        g_fn = permute(g_fn,[2 1 3]);
        
        if ~isempty(varargin)
            
            if numel(varargin) == 2
                rmzero = 1;
                uniqueCombo = 1;
            else
                if strcmp(varargin{1},'unique')
                    uniqueCombo = 1;
                    rmzero = 0;
                else
                    rmzero =1;
                    uniqueCombo = 0;
                end
            end

            if rmzero
            % remove all combinations which contain a zero vector
            tmp = reshape(g_fn,l,[]);
            zero_multis = sum(tmp,1)== 0;
            zero_combos = find(sum(reshape(zero_multis,n,[]),1));
            nonzero_combos = sum(reshape(zero_multis,n,[]),1) == 0;
            %'trivial'
            %g_fn(:,:,zero_combos)
            
            g_fn = g_fn(:,:,nonzero_combos);
            end
            
            if uniqueCombo   % return only unique combinations
                
                % sort pages
                g_fn = g_fn;
                %'l in first dimension'
                for i = 1:size(g_fn,3)
                    combo = sortrows(g_fn(:,:,i).');
                    g_fn(:,:,i) = combo.';
                end
                
                
                % extract unique pages
                [n,m,p]=size(g_fn);
                g_fn=reshape(g_fn,n,[],1);
                g_fn=reshape(g_fn(:),n*m,[])';
                
                [g_fn,~,page_idx] = unique(g_fn,'rows','stable');
                
                g_fn = g_fn.';
                
                % count occurences of the permutations
                n_occurrences = accumarray(page_idx, 1);
                
                g_fn = reshape(g_fn,n,m,[]);
                
                multiplicity{f,n_run} = n_occurrences;
            end
        end
        g{f,n_run}  = g_fn;
        g_comb(f,n_run) = size(g_fn,3);
    end
    
    n_run = n_run +1;
end

if ~isempty(varargin) && uniqueCombo
    varargout{1} = multiplicity;
end
end