function dfnl = compute_dfnldx(obj,x,xd)
% COMPUTE_DFNLDX This function computes the Jacobian of  nonlinear internal 
% force with respect to the displacements x in a second-order
% mechanical system. Currently, we do not treat velocity dependent
% nonlinearity.

assert(obj.order == 2, ' dfnldx can only be computed for second-order systems')

dfnl = sparse(obj.n,obj.n);
switch obj.Options.Intrusion
    case 'full'
        for j = 1:length(obj.fnl)
            if isfield(obj.fnl(j),'subs')
                dfnl = dfnl + expand_tensor_derivative(obj.fnl(j),x);
            else
                dfnl = dfnl + expand_multiindex_derivative(obj.fnl(j),x);
            end
        end

    case 'semi'
        error('implementation for semi-intrusive computation is not available')

    case 'none'
        dfnl = obj.dfnl_non(x);

    otherwise
        error('options for compute fnl are none, semi or full')
end