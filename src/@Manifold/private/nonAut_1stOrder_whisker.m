function [W1,R1,varargout] = nonAut_1stOrder_whisker(obj,W0,R0,input_data)
% NONAUT_1STORDER_WHISKER This function computes the non-autonomous SSM using   
% 
% the first order system algorithm for SSM computation.
%
% [W1,R1] = NONAUT_1STORDER_WHISKER(obj,W0,R0,input_data)
%
% obj:        SSM class object
% W0:         autonomous SSM coefficients
% R0:         autonomous RD coefficients
% input_data: data struct containing necessary information for computation
%
% W1:         non-autonomous SSM coefficients
% R1:         non-autonomous RD coefficients
%
% See also: NONAUT_2NDTORDER_WHISKER, COMPUTE_PERTURBER_WHISKER

% Unpack variables
order  = input_data.order;
nKappa = input_data.nKappa;
N      = input_data.N; 
l      = input_data.l;   


% Structs for storing coefficients,
Fext   = obj.System.Fext;

if obj.Options.contribNonAuto % whether to ignore higher order
    [W1,R1,input_data.kappas,input_data.Fext_ord] = nonAut_struct_setup(l,N,nKappa,order,Fext);
else %only zeroth order autonomous coefficients
    [W1,R1,input_data.kappas,input_data.Fext_ord] = nonAut_struct_setup(l,N,nKappa,0,Fext);
end


% Variables needed for zeroth order computation
data         = input_data;
data.NonAuto = obj.Options.contribNonAuto;
data.solver  = obj.Options.solver;

% System variables
data.A      = obj.System.A;           % A matrix
data.B      = obj.System.B;           % B matrix
data.FEXT   = Fext;

% Manifold variables
data.W_M    = obj.E.adjointBasis ;    % Right eigenvectors of the modal subspace
data.V_M    = obj.E.basis;            % Left eigenvectors of the modal subspace
data.Lambda_M_vector = obj.E.spectrum;


% #########################################################################
%          leading order contributions 
[W1,R1,soltime,nRHS] = nonAut_1stOrder_leadTerms(W1,R1,data);
varargout{1} = nRHS;
if order > 0 || ~obj.Options.contribNonAuto
    % save information about runtime
    obj.solInfoNonAut.timeEstimate(1) = obj.solInfoNonAut.timeEstimate(1) + soltime;
end
% #########################################################################
%          higher order contributions 
if order >0 && obj.Options.contribNonAuto 
    [W1,R1] = nonAut_1stOrder_highTerms(obj,W0,R0,W1,R1,input_data);
end
end

