function [cell_revlex] = coeffs_lex2revlex(cell_lex,string)
% COEFFS_LEX2REVLEX 
%
% Converts input coefficients from lexicographic ordering to reverse
% lexicographic ordering (and also reverse) 
%
% [cell_revlex] = COEFFS_LEX2REVLEX(cell_lex,string)
%
% cell:     cell containing coefficients in (reverse) lexicographic ordering
% string:   specifies which type of coefficients are to be converted
%
% cell_revlex:
%           cell with coefficients in opposite order to the input ordering
%
% See also: COEFFS_CONJ2FULL, COEFFS_CONJ2LEX, COEFFS_SETUP
%%

%%
order = numel(cell_lex);

switch string
    case 'CompCoeff'

        cell_revlex = cellfun(@(x) flip(x,2), cell_lex, 'UniformOutput',false);
        
    case 'TaylorCoeff'
        if isfield(cell_lex,'ind') % For use in nonaut computation
            cell_revlex = repmat(struct('coeffs',[],'ind',[]),order,1);
            for i = 1:order
                cell_revlex(i).coeffs = flip(cell_lex(i).coeffs,2);
                cell_revlex(i).ind    = flip(cell_lex(i).ind).';
            end
            
        else % For use in aut computation
            cell_revlex = repmat(struct('coeffs',[]),order,1);
            for i = 1:order
                cell_revlex(i).coeffs = flip(cell_lex(i).coeffs,2);
            end
        end
end