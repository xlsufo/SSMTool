function dF = set_dF(obj,dF)
% sets nonlinearity in second order form in tensor format
if ~iscell(dF) % multi-index input
    dF = set_F_tensor(dF);
end
end