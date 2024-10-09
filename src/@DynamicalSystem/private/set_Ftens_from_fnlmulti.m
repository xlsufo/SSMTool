function [F] = set_Ftens_from_fnlmulti(obj)
% Converts second order multi-index nonlinearity to first order tensor
% nonlinearity

d = length(obj.fnl) + 1;
F = cell(1,d);
F{1} = sptensor(obj.A);

for j = 2:d
    sizej = obj.N*ones(1,j+1);
    if isempty(obj.fnl(j-1)) || isempty(obj.fnl(j-1).coeffs)
        F{j} = sptensor(sizej);
    else
        [fnl_t] = multi_index_to_tensor(obj.fnl(j-1).coeffs,obj.fnl(j-1).ind);
        subsj = fnl_t.subs;
        valsj = -fnl_t.vals;
        if obj.order==1
            valsj = -valsj;
        end
        F{j} = sptensor(subsj,valsj,sizej);
    end
    
end
end