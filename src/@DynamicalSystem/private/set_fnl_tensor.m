function [fnl_tensor] = set_fnl_tensor(obj,fnl)
% Sets second order nonlinear force in tensor format

d = length(fnl) + 1;
fnl_tensor = cell(1,d);

for j = 2:d
    if size(fnl(j-1).ind,2) == obj.N
        sizej = [obj.n, obj.N*ones(1,j)];
    else
        sizej = obj.n*ones(1,j+1);
    end
    if isempty(fnl(j-1)) || isempty(fnl(j-1).coeffs)
        fnl_tensor{j-1} = sptensor(sizej);
    else
        [fnl_t] = multi_index_to_tensor(fnl(j-1).coeffs,fnl(j-1).ind);
        subsj = fnl_t.subs;
        valsj = fnl_t.vals;
        if obj.order==1
            valsj = -valsj;
        end
        fnl_tensor{j-1} = sptensor(subsj,valsj,sizej);
    end
    
end
end