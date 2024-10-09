function [RES,I_minu,I_subt] = multi_subtraction(MINUEND,SUBTRAHEND,type)
% MULTI_SUBTRACTION  This function subtracts multi-indices.
%
% The output is a matrix, containing all columns of $\texttt{MINUEND}$ subtracted
% by all columns of $\texttt{SUBTRAHEND}$ that do not lead to columns with negative entries.
% The matrix $\texttt{I\_subt}$ contains in position $j$ the column index of $\texttt{SUBTRAHEND}$
% corresponding to the column $j$ of $\texttt{RES}$. The same holds for $\texttt{I_minu}$.
%
% Example:
%
% $\texttt{I\_subt(3) = 1}$, $\texttt{I\_minu(3) = 5}$ implies that  $\texttt{RES(:,3)}
% = \texttt{MINUEND(:,5)} - \texttt{SUBTRAHEND(:,1)}$. 
%
% [RES,I_minu,I_subt] = MULTI_SUBTRACTION(MINUEND,SUBTRAHEND,type)
%
% MINUEND:    matrix containing the minuend multi-indices in its columns
% SUBTRAHEND: matrix containing the subtraend multi-indices in its columns
% type:       'Arbitrary' - generic multi-indices that are subtracted
%             'Identity'  - Subtrahend is the identity matrix  
%             'Unit'      - Subtrahend is collection of unit vectors
%
% RES:        matrix containing all subtracted multi-index pairs with
%             nonnegative entries
% I_minu:     indices minuends leading to RES
% I_subt:     indices subtrahends leading to RES
%
% See also: MULTI_ADDITION, MULIT_NSUMK, MULTI_INDEX_2_ORDERING



%%
% different functionalities for parametrisation multi indices and forcing multi
% indices
switch type
    case 'Arbitrary'
        minu_sz = size(MINUEND,2);
        subt_sz = size(SUBTRAHEND,2);
        
        % Index of the subtrahends
        I_subt = reshape(repmat(1:subt_sz,minu_sz,1),1,[]);
        % Index of the minuends
        I_minu = repmat(1:minu_sz,1,subt_sz);
        RES = MINUEND(:,I_minu)-SUBTRAHEND(:,I_subt);
        
        
        [~,neg]    = find(RES<0);
        RES(:,neg)   = []; %delete columns with negative entries
        I_subt(:,neg) = [];
        I_minu(:,neg) = [];
        
        
    case 'Identity'
        %subtrahend is always identity matrix

        [I_minu,I_subt] = find(MINUEND.');
        RES =  (MINUEND(:,I_minu)-SUBTRAHEND(:,I_subt));
        
    case 'Unit' 
        % subtrahend is collection of unit vectors
        [I_subt_col,I_subt_row] = find(SUBTRAHEND.'); % contains index of unit vectors
        [I_minu,I_subt] = find(MINUEND(I_subt_row,:).'); % finds nonzero entries in Minuend in relevant rows
        
        RES = MINUEND(:,I_minu) - SUBTRAHEND(:,I_subt_col(I_subt));
        I_subt = I_subt_col(I_subt);
end

% Output indices as row vector
if size(I_minu,1)>1
    I_minu = I_minu.';
end

if size(I_subt,1)>1
    I_subt = I_subt.';
end

end