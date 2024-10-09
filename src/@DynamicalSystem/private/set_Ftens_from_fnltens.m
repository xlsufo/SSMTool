function [F] = set_Ftens_from_fnltens(obj)
% Converts second order multi-index nonlinearity to first order tensor
% nonlinearity

d = length(obj.fnl)+1;
F = cell(1,d);


for j = 2:d
    sizej = obj.N*ones(1,j+1);
    if isempty(obj.fnl{j-1})
        F{j} = sptensor(sizej);
    else
        subsj = obj.fnl{j-1}.subs;
        valsj = -obj.fnl{j-1}.vals;
        if obj.order==1
            valsj = -valsj;
        end
        F{j} = sptensor(subsj,valsj,sizej);
    end
    
end
end