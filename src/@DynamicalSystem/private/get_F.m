function F = get_F(obj,obj_F)

if isempty(obj_F) && ~isempty(obj.fnl)
    F = get_F_from_fnl(obj);
    
    obj.F = F;
else
    F = obj_F;
end
end