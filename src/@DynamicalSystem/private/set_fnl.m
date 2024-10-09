function fnl = set_fnl(obj,fnl)

% sets nonlinearity in second order form in tensor format
if ~iscell(fnl) % multi-index input
    fnl = set_fnl_tensor(obj,fnl);
end

end