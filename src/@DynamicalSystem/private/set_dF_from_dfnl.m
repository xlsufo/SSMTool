function dF = set_dF_from_dfnl(obj)
% Converts second order multi-index nonlinearity to first order tensor
% nonlinearity

d = length(obj.dfnl)+1;
dF = cell(1,d);

for j = 2:d
    sizej = obj.N*ones(1,j+1);
    if isempty(obj.dfnl{j-1})
        dF{j} = sptensor(sizej);
    else
        subsj = obj.dfnl{j-1}.subs;
        valsj = -obj.dfnl{j-1}.vals;
        if obj.order==1
            valsj = -valsj;
        end
        dF{j} = sptensor(subsj,valsj,sizej);
    end
    
end
end