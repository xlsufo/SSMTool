function fnl = compute_fnl(obj,x,xd)
% COMPUTE_FNL We compute the nonlinear internal force in a second-order
% mechanical system. 

assert(obj.order == 2, ' fnl can only be computed for second-order systems')

fnl = zeros(obj.n,1);
if isempty(obj.fnl) && strcmp(obj.Options.Intrusion,'none')
    if obj.nl_input_dim == obj.N
        fnl = obj.fnl_non([x;xd]);
    else
        fnl = obj.fnl_non(x);
    end
else
    for j = 1:length(obj.fnl)
        if obj.nl_input_dim == obj.N % check if the nonlinearity is velocity dependent as well
            z = [x;xd];
        else
            z = x;
        end
        if iscell(obj.fnl)
            % expand_tensor does not support vectorization
            % Fnl = Fnl + expand_tensor(obj.F{j},z);
            fj  = tensor_to_multi_index(obj.fnl{j});
            fnl = fnl + expand_multiindex(fj,z);
        else
            fnl = fnl + expand_multiindex(obj.fnl(j),z);
        end
    end

end