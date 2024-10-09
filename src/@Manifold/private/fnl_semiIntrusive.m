function [Fk] = fnl_semiIntrusive(F,W,nl_order,K,data)
% FNL_SEMIINTRUSIVE This function computes the nonlinear contributions to the invariance
%
% equation for multi indices in set K, given a semi-intrusive function handle. 
% This function assumes the displacement parametrisation to be in the first
% n entries of the full SSM parametrisation, ie a first order system of the
% form z = [x ; dot(x)]
%
% [Fk]= FNL_SEMIINTRUSIVE(F,W,nl_order,K,data)
%
% F:        semi-function handle for the nonlinearity
% W:        autonomous SSM coefficients
% nl_order: order of nonlinearity considered
% K:        set of multi indices at current order of computation
% data:     struct containing necessary information for the computation
%
% Fk:       F composed with aut SSM parametrisation W, contributions at
%           order K
%
% See also: COHOMOLOGICAL_SOLUTION, FNL_NONINTRUSIVE

switch data.ordering
    case 'revlex'
        [Fk] = nl_revlex(F,W,nl_order,K,data);
    case 'conjugate'
        [Fk] = nl_conj(F,W,nl_order,K,data);
end
end

function [Fk] = nl_conj(F,W,nl_order,K,multi_input)
N = size(W(1).coeffs,1);
z_k = size(K,2);
Fk = sparse(N,z_k);

l = multi_input.l;
ordering = multi_input.ordering;
revlex2conj = multi_input.revlex2conj;
% Conjugate center indices
Z_cci = multi_input.Z_cci;

% dimensionality of input vector for nonlinearity
mx_idx = multi_input.nl_input_dim;

for m = 1:z_k

    % Redundantly also computes combos with zero vectors
    % Assumes symmetric function handles wrt. input vectors
    switch multi_input.symmetry
        case true
            [g,~,~,~,perms] = multi_nsumk(nl_order,K(:,m),'unique'); 
            g = g{1};
            perms = perms{1};
        case false
            [g,~,~,~] = multi_nsumk(nl_order,K(:,m)); 
             g = g{1};
    end
    
    Fk_m = sparse(N,1);
    
    % Loop over all partitions of m
    for i = 1:size(g,3)
        h_abs = sum(g(:,:,i));
        
        % Check if any multi-index in set is all zero
        if any(h_abs == 0)
            continue
        else
            
            h_idx = multi_index_2_ordering(g(:,:,i),ordering,revlex2conj);
            
            % Check for all zero SSM coefficient vector
            vectors = cell(nl_order,1);
            emptyflag = false;
            %[vectors{:}] = deal(W(h_abs).coeffs);
            for j = 1:nl_order
                
                % check if conjugate multi-index is to be used                
                if h_idx(j) > Z_cci(h_abs(j))
                    
                    z_ord = nchoosek(h_abs(j)+l-1,l-1);
                    
                    vectors{j} = conj(W(h_abs(j)).coeffs(1:mx_idx,z_ord -h_idx(j) +1 ));
                    
                else
                    vectors{j} = W(h_abs(j)).coeffs(1:mx_idx,h_idx(j));
                end
                
                if isempty(vectors{j})
                    emptyflag = true;
                    continue
                end
                    
            end

            if emptyflag % no contribution of nl
                continue
            end
                                         
           
            switch multi_input.symmetry
                case true
                    Fk_m = Fk_m+ perms(i)*F( vectors);
                    
                case false                  
                    Fk_m = Fk_m+  F( vectors);                    
            end
        end
    end

    Fk(:,m) = Fk_m;
end

Fk = double(Fk);

end

function [G] = nl_revlex(F,W,nl_order,K,multi_input)
N = size(W(1).coeffs,1);
z_k = size(K,2);
G = sparse(N,z_k);

ordering = multi_input.ordering;

% dimensionality of input vector for nonlinearity
mx_idx = multi_input.nl_input_dim;

for m = 1:z_k

    % Redundantly also computes combos with zero vectors
    switch multi_input.symmetry
        case true
            [g,~,~,~,perms] = multi_nsumk(nl_order,K(:,m),'unique'); 
            g = g{1};
            perms = perms{1};
        case false
            [g,~,~,~] = multi_nsumk(nl_order,K(:,m)); 
             g = g{1};
    end
    
    Gm = sparse(N,1);
    % Loop over all partitions of m
    for i = 1:size(g,3)
        h_abs = sum(g(:,:,i));
        
        % Take out combinations which have a zero vector
        if any(h_abs == 0)
            continue
        else
            
            % get indices of multi-index set
            h_idx = multi_index_2_ordering(g(:,:,i),ordering,[]);
            
            vectors = cell(nl_order,1);
            %[vectors{:}] = deal(W(h_abs).coeffs);
            
            % Check for all zero SSM coefficient vector
            emptyflag = false;
            
            for j = 1:nl_order          

                    vectors{j} = W(h_abs(j)).coeffs(1:mx_idx,h_idx(j));
                    
                    if isempty(vectors{j})
                        emptyflag = true;
                        continue
                    end
            end
            
            if emptyflag % no contribution of nl
                continue
            end
            
            switch multi_input.symmetry
                case true
                    Gm = Gm + perms(i)*F( vectors);
                case false
                    Gm = Gm +F( vectors);
                    
            end
        end
    end

    G(:,m) = Gm;
end

G = double(G);

end