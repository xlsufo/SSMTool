function F_non  = get_F_non(obj,obj_F_non)
% F_NON
% Get method


if obj.order ==1  && ~isempty(obj_F_non)
    F_non = obj_F_non;

elseif (obj.order == 2 && ~isempty(obj.fnl_non)) || (isempty(obj_F_non) && ~isempty(obj.fnl_non)) % second order system, or second order and computation carried out first order
    F_non = set_F_non(obj);

else
    F_non = [];
end
end