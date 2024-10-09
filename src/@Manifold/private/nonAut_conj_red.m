function [redConj,mapConj] = nonAut_conj_red(kappa_set,F_kappa)
% NONAUT_CONJ_RED This function detects complex conjugate relations between forcing. 

% For instance, when kappa_set = [1,-1,2,3,-3] and F_kappa = [1;1;2;3;4], 
% it will return redConj = [1,3,4,5] with mapConj = {[1 2],3,4,5}  
%
% [redConj,mapConj] = NONAUT_CONJ_RED(kappa_set,F_kappa)
%
% kappa_set: set of harmonics
% F_kappa:   forcing coefficients
%
% redConj:   Index of non conjugate harmonics in kappa_set
% mapConj:   Index of coefficients in F_kappa which belong to
%            kappa_set(redConj)
%
% See also: NONAUT_1STORDER_LEADTERMS, NONAUT_2NDORDER_LEADTERMS

redConj = [];
mapConj = [];
assert(numel(kappa_set)==numel(unique(kappa_set)),'there exist redundancy in kappa of external forcing');
kappa = kappa_set;
while ~isempty(kappa)
    ka = kappa(1);
    ka_redConj = find(kappa_set==ka);
    redConj = [redConj;ka_redConj];
    % find the conjugate one if it exists
    ka_conj = find(kappa_set==-ka);
    if ~isempty(ka_conj) && norm(conj(F_kappa(:,ka_redConj))-F_kappa(:,ka_conj))<1e-6*norm(F_kappa(:,ka_conj))
        mapConj = [mapConj, {[ka_redConj,ka_conj]}];
        kappa = setdiff(kappa,[ka,-ka],'stable');
    else
        mapConj = [mapConj, {ka_redConj}];
        kappa = setdiff(kappa,ka,'stable');
    end
end
end