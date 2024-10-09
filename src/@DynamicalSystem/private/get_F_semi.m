function F_semi    = get_F_semi(obj,obj_F_semi)
if obj.order ==1 && ~isempty(obj_F_semi)
    F_semi = obj_F_semi;

elseif obj.order == 2 || (isempty(obj_F_semi) && ~isempty(obj.fnl_semi)) % second order system, or second order and computation carried out first order
    F_semi = construct_F_semi(obj);

end
end

function [F_semi]    = construct_F_semi(obj)
% Checks input dimensions of vectors and unifies to f_func
F_semi = cell(numel(obj.fnl_semi)+1,1 );

for j = 2:numel(obj.fnl_semi)+1
    %invecs is a N times nl_order dimensional array containing a input
    %vector in each column
    if ~isempty(obj.fnl_semi{j-1})
    F_semi{j} = @(invecs) [-double(obj.fnl_semi{j-1}(invecs)); sparse(obj.n,1)];
    end
end
end