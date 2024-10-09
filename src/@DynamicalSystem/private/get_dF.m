function dF = get_dF(obj,obj_dF)

if isempty(obj_dF)
    if ~isempty(obj.F)
        dF = construct_dF(obj);
    elseif ~isempty(obj.dfnl)
        dF = set_dF_from_dfnl(obj);
    else
        dF = [];
    end
    set(obj,'dF',dF)
else
    if ~iscell(obj_dF)
        error('Jacobian must be input in tensor format')
    end
    dF = obj_dF;
end

end


function dF = construct_dF(obj)
% Compute jacobians from internal forces, if they are not provided
% explicitly.

dF = {};

for i = 2:numel(obj.F)
    if ~isempty(obj.F(i)) && nnz(obj.F{i})>0
        dim = 2:(i+1);
        dFi = sptensor(size(obj.F{i}));
        for perm_dim = 2:i+1
            dFi = dFi + permute(obj.F{i}, [1,perm_dim, dim(dim~=perm_dim)]);
        end
        dF{i} = dFi;
    end
end

end