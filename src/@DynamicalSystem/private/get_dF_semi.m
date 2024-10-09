function dF_semi   = get_dF_semi(obj,obj_dF_semi)
if obj.order ==1 && ~isempty(obj_dF_semi)
    dF_semi = obj_dF_semi;

elseif obj.order == 2 || (isempty(obj_dF_semi) && ~isempty(objdf.nl_semi)) % second order system, or second order and computation carried out first order
    if isempty(obj_dF_semi)
        dF_semi = construct_dF_semi(obj);
        obj.dF_semi = dF_semi;
    else
        dF_semi = obj_dF_semi;
    end
end
end

function [dF_semi]   = construct_dF_semi(obj)
% Checks input dimensions of vectors and unifies to f_func
dF_semi = cell(numel(obj.dfnl_semi)+1,1 );

for j = 2:numel(obj.dfnl_semi)+1
    %invecs is a N times nl_order dimensional array containing a input
    %vector in each column
    if ~isempty(obj.dfnl_semi{j-1})

            dF_semi{j} = @(invecs) [-double(obj.dfnl_semi{j-1}(invecs)); sparse(obj.n,1)];
    end
end
end