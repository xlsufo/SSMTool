function [R_0i,W_0i,varargout] = Aut_2ndOrder_SSM(WR,Fn,data,Mass,Damp,Stiff)
% AUT_2NDORDER_SSM This function computes the autonomous SSM at order k using   
% 
% the second order system algorithm for SSM computation.
%
% [R_0i,W_0i,varargout] = AUT_2NDORDER_SSM(WR,Fn,data,Mass,Damp,Stiff)
%
% WR:       product of SSM and RD coefficients that contribute at order k
% Fn:       internal force composed with SSM evaluated at order k
% data:     data struct containing necessary information for computation
% Mass:     Mass matrix of system
% Damp:     Damping matrix of system
% Stiff:    Stiffness matrix of system
%
% R_0i:     order k autonomous reduced dynamics
% W_0i:     order k autonomous SSM coefficients
% varargout:information concerning the computation times.
%
% See also: AUT_1STORDER_SSM, COHOMOLOGICAL_SOLUTION

%retrieve inputs
[z_k,l,N,K,Lambda_M_vector,solver,reltol] = deal(data.z_k,data.l,data.N,data.K,data.Lambda_M_vector,data.solver,data.reltol);
THETA  = data.W_M(1:(N/2),:);
PHI    = data.V_M(1:(N/2),:);


Ym     = -( Damp * WR(1:(N/2),:) + Mass * WR((N/2+1):end,:)) - Fn(1:(N/2),:);
Vm     = - WR(1:(N/2),:);


Lambda_K         = sum(K.*Lambda_M_vector);

[I,F] = Aut_resonant_terms(Lambda_M_vector,Lambda_K,reltol); % F contains multi-index positions

%% Analytic Reduced dynamics
% {
[Rk] = Aut_2ndOrder_RedDyn(I,F,THETA,PHI,Damp,Lambda_K,Lambda_M_vector, Mass,Vm,Ym,l,z_k);

w_0i    = zeros(N/2,z_k);
w_0idot = zeros(N/2,z_k);


eqtoc = 0;
for f = 1:z_k
    L_k = ( Mass * ((Lambda_K(f) + Lambda_M_vector.') .* PHI) + Damp*PHI ) * Rk(:,f);
    L_k = L_k + Lambda_K(f)*Mass*Vm(:,f) + Ym(:,f);
        
    C_k = -(Stiff + Lambda_K(f)*Damp + Lambda_K(f)^2 *Mass );
    eqf = tic;
    if any(L_k)
        w_0i(:,f)    = solveinveq(C_k,L_k,solver);
    end
    eqtoc = eqtoc+toc(eqf);
    
    w_0idot(:,f) = Lambda_K(f) * w_0i(:,f) + PHI * Rk(:,f) + Vm(:,f);

end

%}
%% Bordered approach
%{
w_0i_test    = zeros(N/2,z_k);
w_0idot_test = zeros(N/2,z_k);
Rk_test = zeros(l,z_k);

eqtoc = 0;
for f = 1:z_k
    
    L_k = Lambda_K(f)*Mass*Vm(:,f) + Ym(:,f);
    
    C_k = -(Stiff + Lambda_K(f)*Damp + Lambda_K(f)^2 *Mass );
    
    eqf = tic;
    if ismember(f,F) && any(L_k)

        THETA_f = THETA(:,I(F==f));
        PHI_f   = PHI(:,I(F==f));
                
        C11 = -C_k;
        C12 = ( Mass * ((Lambda_K(f) + Lambda_M_vector(I(F==f)).') .* PHI_f) + Damp*PHI_f );
        C21 = THETA_f' .* ((Lambda_K(f) + Lambda_M_vector(I(F==f))) )*Mass + THETA_f' * Damp ;
        C22 =  THETA_f' *  Mass * PHI_f ;
        
        C_k = [ C11 , C12 ; C21 , C22 ];
        Lend = - THETA_f' * Mass * Vm(:,f);
        L_k = [-L_k; Lend];
        

        sol  =  solveinveq(C_k,L_k,solver);
        w_0i_test(:,f) = sol(1:end-1);
        Rk_test(I(F==f),f) = sol(end);
    
    elseif any(L_k)
        w_0i_test(:,f)    = solveinveq(C_k,L_k,solver);
    end
    
    eqtoc = eqtoc+toc(eqf);
    
    w_0idot_test(:,f) = Lambda_K(f) * w_0i_test(:,f) + PHI * Rk(:,f) + Vm(:,f);

end
%}

%% Prepare output
W_0i       = [w_0i;w_0idot];
R_0i       = Rk;
rdtoc = 0;


varargout{1} = eqtoc;
varargout{2} = rdtoc;
end