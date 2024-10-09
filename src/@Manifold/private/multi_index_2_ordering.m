function [I] = multi_index_2_ordering(m,ordering,lex2conj)
% MULTI_INDEX_2_ORDERING This function computes the position of a multi
%
% index in an ordered set of all multi-indices of the same order.
% Given a multi- index array containing multi-indices in its columns this function computes
% the position of the multi-indices in a ordered set of all multi-index
% vectors that have the same order as each of the multi-indices. The ordering can be
% 'lex' for lexicographical or 'antilex for anti-lexicographical or 'conj'
% for conjugate ordering, documentation on conjugate ordering in the file
% $\textit{Explicit Kernel
% Extraction and Proof ofSymmetries of SSM Coefficients - Multi-Indexversion}$
% which can be found in the Documentation branch.
%
% [I] = MULTI_INDEX_2_ORDERING(m,ordering,lex2conj)
%
% m:        multi-inidces
% ordering: string that determines which order the set has
% lex2conj: index set to convert from lexicographic to conjugate ordering,
%           optional parameter, only needed for conjugate ordering.
%
% I:        Set containing all the indices of multi-indices in m
%
% See also:  MULTI_ADDITION, MULIT_NSUMK, MULTI_SUBTRACTION


%%
if isempty(m)
    I = [];
elseif size(m,1) == 1 % 1D SSM
    sz   = size(m,2);
    I    = ones(1,sz);    
else
    l    = size(m,1);
    sz   = size(m,2);
    mabs = sum(m,1);
    I    = ones(1,sz);
    for ij  = 1:sz
        bij = 0;
        if m(1,ij) ~= 0
            for i = 0:m(1,ij)-1
                bij = bij+ nchoosek(mabs(ij)-i+l-1-1,l-1-1);
            end
        end
        
        for j = 2:(l-1)
            if m(j,ij) ~= 0
                for i  = 0:(m(j,ij)-1)
                    bij = bij+ nchoosek(mabs(ij)-sum(m(1:j-1,ij)) ...
                        -i+l-j-1,l-j-1);
                end
            end
        end
        I(ij) = bij+1;
    end
    
    %if order is reverse-lexicographical, then the position is
    % the mirrored value about the centre of the array.
    switch ordering
        case 'lex'
            %Ordering is corrects
        
        case 'revlex'
            if all(mabs == mabs(1))
                I = nchoosek(mabs(1)+l-1,l-1)+1-I;
            else
                for i = 1:sz
                    I(i) = nchoosek(mabs(i)+l-1,l-1)+1-I(i);
                end
            end
     
    % if order is conjugate, convert to that         
        case 'conjugate'
            if all(mabs == mabs(1))
                I = nchoosek(mabs(1)+l-1,l-1)+1-I;    
            else
                for i = 1:sz
                    I(i) = nchoosek(mabs(i)+l-1,l-1)+1-I(i);
                    
                end
            end

            for i = 1:sz
                idx = lex2conj{mabs(i)};
                I(i)= find(idx==I(i));
            end
    end
end
end

