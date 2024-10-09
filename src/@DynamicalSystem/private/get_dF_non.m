function dF_non  = get_dF_non(obj,obj_dF_non)
% GET_DF_NON
% Get method
if obj.order ==1 && ~isempty(obj_dF_non)
    dF_non = obj_dF_non;

elseif obj.order == 2 || (isempty(obj_dF_non) && ~isempty(obj.dfnl_non)) % second order system, or second order and computation carried out first order
    dF_non = set_dF_non(obj);
    set(obj,'dF_non',dF_non)
end

end