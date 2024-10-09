function [FGs] = nonAut_Fext_plus_Fnl(NL,FEXT,data,k,W0,W1,nKappas)
% NONAUT_FEXT_PLUS_FNL  This function computes the contribution to higher order 

% non-autonomous cohomological equations of the nonlinear internal forces
% and the possibly linear or nonlinear external excitation.
%
% [FGs] = NONAUT_FEXT_PLUS_FNL(NL,FEXT,data,k,W0,W1,nKappas)
%
% NL:       information about the internal forces, either in tensor form or
%           as function handle
% FEXT:     struct array containing the external forcing
% data:     data struct containing necessary information for computation
% k:        current order of SSM computation
% W0:       autonomous SSM coefficients
% W1:       non-autonomous SSM coefficients
% nKappas:  number of harmonics in forcing
%
% FGs:      internal and external forces composed with SSM coefficients for
%           all harmonics.
%
% See also: NONAUT_1STORDER_HIGHTERMS, NONAUT_2NDORDER_HIGHTERMS

z_k   = nchoosek(k+data.l-1,data.l-1);
FGs   = repmat(struct('val' ,[]),nKappas, 1);
Force = repmat(struct('val' ,sparse(data.N,z_k)),nKappas, 1);
DF    = repmat(struct('val' ,sparse(data.N,z_k)),nKappas, 1);

if data.l > 1
    K   = flip(sortrows(nsumk(data.l,k,'nonnegative')).',2); %order k multi-indices
else
    K = k;
end
%% Assemble linear & nonlinear terms
for n = 2:k+1
    % FORCING
    %sum to k+1 since index starts at 0 for k=0
    for i = 1:nKappas
    if  n <= data.Fext_ord(i) && ~isempty(FEXT.data(i).F_n_k(n).coeffs)
        F_coeff = FEXT.data(i).F_n_k(n).coeffs;
        F_ind   = FEXT.data(i).F_n_k(n).ind.';
        Force(i).val   = Force(i).val + F_coeff * fnl_intrusive(F_ind,K, data); %fnl_intrusive computes composition
        
    end
    end
    
    switch data.mode
        case 'SemiIntrusiveF'
            if n <= data.nl_ord && ~isempty(NL{n})
                DFi = dfnl_semiIntrusive(NL{n},n,W0,W1,K,data,nKappas);
                
                for i = 1:nKappas
                    DF(i).val  = DF(i).val + DFi(i).val  ;
                end
            end
            
        case 'IntrusiveF'
            if n <= data.nl_ord && ~isempty(NL{n})
                DFi = dfnl_semiIntrusive(NL{n},n,W0,W1,K,data,nKappas);
                
                for i = 1:nKappas
                    DF(i).val  = DF(i).val + DFi(i).val  ;
                end
            end
            
        case 'NonIntrusiveF'
            if n <= data.nl_ord && ~isempty(NL) 
                DFi = dfnl_nonIntrusive(NL,n,W0,W1,K,data,nKappas);                                              
                
                for i = 1:nKappas
                    DF(i).val  = DF(i).val + DFi(i).val  ;
                end
            end            
    end
end


for i = 1:nKappas
    FGs(i).val = Force(i).val + DF(i).val;
end

end
