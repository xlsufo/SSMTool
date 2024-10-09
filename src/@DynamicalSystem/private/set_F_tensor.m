function [F_tensor] = set_F_tensor(obj)
% Sets second order nonlinear force in tensor format

d = length(obj.F);
F_tensor = cell(1,d);

for j = 2:d
    sizej = obj.N*ones(1,j+1);
    if isempty(obj.F(j)) || isempty(obj.F(j).coeffs)
        F_tensor{j} = sptensor(sizej);
    else
        [F_t] = multi_index_to_tensor(obj.F(j).coeffs,obj.F(j).ind);
        subsj = F_t.subs;
        valsj = F_t.vals;
        if obj.order==1
            valsj = -valsj;
        end
        F_tensor{j} = sptensor(subsj,valsj,sizej);
    end
    
end
end