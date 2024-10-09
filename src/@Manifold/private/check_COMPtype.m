function [COMPtype, DStype] = check_COMPtype(obj)
% CHECK_COMPTYPE
%
% Check if the SSM computation should be carried out using the first order
% formalism. This may be due to user input, a complex valued dynamical
% system that is underlying, for 1-dimensional SSMs.
% Furthermore this function checks whether the symmetries in conjugate
% multi-indices can be taken advantage of, which is only possible for real
% dynamical systems
%
% By default, COMPtype is set to 'first'. In general, for second order
% systems the user may choose, which computation approach he desires to
% use. Some cases do however demand the first order computation. In the
% following we ensure, that in these cases no second order computations are
% performed.
%
% [COMPtype, DStype] = CHECK_COMPTYPE(obj)
%
% obj:      instance of the DS class
%
% COMPtype: string that determines which algorithm is used
% DStype:   string that determines whether DS is real or complex/SSM is 1D
%
% See also: DynamicalSystem


DStype = checkDStype(obj);

if obj.System.order == 1
    % First order system -> first order computation
    COMPtype  = 'first'; 
    
elseif strcmp(obj.Options.COMPtype,'second')
    % Second order system with second order computation.
    
    if strcmp(DStype,'complex')
        % Complex dynamical systems are only supported in first order
        % computation
        
        fprintf('\n Second order SSM computation only supported for real systems, using first order algorithm. \n')
        COMPtype = 'first';
        
    else
        COMPtype = 'second'; %second order computation can be performed
        
    end
    
else
    % Second order system with default first order computation
    COMPtype  = 'first';
end


if obj.dimManifold == 1
    % System is not complex but computational routine works exactly the
    % same for 1D SSM as for complex DS
    fprintf('\n Computation of 1D SSM is performed using first order algorithm. \n')

    COMPtype  = 'first';
    DStype    = 'complex';
end

obj.Options.COMPtype = COMPtype;
obj.System.Options.DStype   = DStype;


end

function [DStype]   = checkDStype(obj)

%% Check if System is complex valued
if ~isreal(obj.System.A) || ~isreal(obj.System.B)
    DStype = 'complex';
    
    fprintf('\n System matrices are complex valued. \n')    
    
%elseif ~checkF(obj)
%    DStype = 'complex';
    
%    fprintf('\n Internal forces are complex valued. \n')
else
    DStype = 'real';
end
end

function [indicator] = checkF(obj)


d = length(obj.System.F) ;
for j = 1:d
    if ~isreal(obj.System.F(j).coeffs)
        indicator = false;
        return
    end
end

indicator = true;
end
