function dfnl = get_dfnl(obj,obj_dfnl)
if ~isempty(obj_dfnl)
    if ~iscell(obj_dfnl)
        error('Jacobian must be input in tensor format')
    end
    dfnl = obj_dfnl;
elseif isempty(obj_dfnl)
    dfnl = construct_dfnl(obj);
    obj.dfnl= dfnl;
else
    dfnl = [];
end

end

function dfnl = construct_dfnl(obj)
% Compute jacobians from internal forces, if they are not provided
% explicitly.

dfnl = {};

for i = 1:numel(obj.fnl)
    if ~isempty(obj.fnl{i})
        dim = 2:(i+2);
        dfnli = sptensor(size(obj.fnl{i}));
        for perm_dim = 2:i+2
            dfnli= dfnli + permute(obj.fnl{i}, [1,perm_dim, dim(dim~=perm_dim)]);
        end
        dfnl{i}=dfnli;
    end
end


end