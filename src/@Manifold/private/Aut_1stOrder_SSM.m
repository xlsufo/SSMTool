function [R_0i,W_0i,varargout] = Aut_1stOrder_SSM(RHS, data, A,B)
% AUT_1STORDER_SSM This function computes the autonomous SSM at order k using   
% 
% the second order system algorithm for SSM computation.
%
% [R_0i,W_0i,varargout] = AUT_2NDORDER_SSM(WR,Fn,data,Mass,Damp,Stiff)
%
% RHS:      right hand side of invariance equation at order k
% data:     data struct containing necessary information for computation
% A:        system matrix
% B:        system matrix
%
% R_0i:     order k autonomous reduced dynamics
% W_0i:     order k autonomous SSM coefficients
% varargout:information concerning the computation times.
%
% See also: AUT_2NDORDER_SSM, COHOMOLOGICAL_SOLUTION

[z_k,l,N, K, Lambda_M_vector,W_M,V_M,reltol,solver] = deal(data.z_k,data.l,data.N, data.K, data.Lambda_M_vector,data.W_M,data.V_M,data.reltol,data.solver);
                    

% Extract the near kernel of the coefficient matrix
% coordinate directions do not change - we use the evals as in rev. lex
% ordering - lambda_i has to be multiplied with i-th entry of a multi-index
K_Lambda         = sum(K.*Lambda_M_vector);

rd = tic;

[R_0i,RHS] = Aut_1stOrder_RedDyn(z_k,Lambda_M_vector,K_Lambda,W_M,reltol,RHS,V_M,B);

rdtoc = toc(rd);

W_0i    = zeros(N,z_k);
RHS    = reshape(RHS,N,z_k);

% Solve the linear system for the SSM-coefficients
eqtoc = 0;

for f = 1:z_k
    if any(RHS(:,f))
        C_k        = B*K_Lambda(f)-A;
        eqf = tic;
        
        W_0i(:,f) = solveinveq(C_k,RHS(:,f),solver);
        eqtoc = eqtoc+toc(eqf);
        
    end
end

R_0i       = reshape(R_0i,l,[]);

varargout{1} = eqtoc;
varargout{2} = rdtoc;
end