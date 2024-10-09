function [W1,R1] = nonAut_2ndOrder_whisker(obj,W0,R0,input_data)
% NONAUT_2NDORDER_WHISKER This function computes the non-autonomous SSM using   
% 
% the second order system algorithm for SSM computation.
%
% [W1,R1] = NONAUT_2NDORDER_WHISKER(obj,W0,R0,input_data)
%
% obj:        SSM class object
% W0:         autonomous SSM coefficients
% R0:         autonomous RD coefficients
% input_data: data struct containing necessary information for computation
%
% W1:         non-autonomous SSM coefficients
% R1:         non-autonomous RD coefficients
%
% See also: NONAUT_1STORDER_WHISKER, COMPUTE_PERTURBER_WHISKER

% Unpack variables
order  = input_data.order;
nKappa = input_data.nKappa;
N      = input_data.N; 
n      = N/2;
l      = input_data.l;   
FEXT   = obj.System.Fext;

% Structs for storing coefficients,
if obj.Options.contribNonAuto % whether to ignore higher order
    [W1,R1,input_data.kappas,input_data.Fext_ord] = nonAut_struct_setup(l,N,nKappa,order,FEXT);
else %only zeroth order autonomous coefficients
    [W1,R1,input_data.kappas,input_data.Fext_ord] = nonAut_struct_setup(l,N,nKappa,0,FEXT);
end

% Variables needed for zeroth order computation
data         = input_data;
data.NonAuto = obj.Options.contribNonAuto;
data.solver  = obj.Options.solver;
% System matrices
data.Mass      = obj.System.M;
data.Damp      = obj.System.C;
data.Stiff     = obj.System.K;
data.FEXT      = FEXT;
% Manifold variables
data.THETA  = obj.E.adjointBasis(1:n,:);
data.PHI    = obj.E.basis(1:n,:);
data.Lambda_M_vector = obj.E.spectrum;

% ########################################################################
%                   leading order contributions
[W1,R1,soltime] = nonAut_2ndOrder_leadTerms(W1,R1,data);

if order > 0 || ~obj.Options.contribNonAuto
    % save information about runtime
    obj.solInfoNonAut.timeEstimate(1) = obj.solInfoNonAut.timeEstimate(1) + soltime;
end

% ########################################################################
%                higher order contributions
if order >0 && obj.Options.contribNonAuto   
   [W1,R1] = nonAut_2ndOrder_highTerms(obj,W0,R0,W1,R1,input_data);    
end

end


