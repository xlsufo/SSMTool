function Fext  = get_Fext(obj,obj_Fext)
% GET_FEXT
% Get method

if obj.order == 1 && ~isempty(obj_Fext)
    Fext = obj_Fext;
elseif obj.order == 1 && isempty(obj_Fext)
    Fext = [];
elseif obj.order == 2 || ( isempty(obj_Fext) && ~isempty(obj.fext)) % second order system with second or first order computation
    Fext.data    = set_Fext(obj);
    Fext.epsilon = obj.fext.epsilon;
end

if isempty (obj_Fext)
    obj.Fext = Fext;
end

end