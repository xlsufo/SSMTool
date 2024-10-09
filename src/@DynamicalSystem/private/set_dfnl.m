function dfnl = set_dfnl(obj,dfnl)
% sets nonlinearity in second order form in tensor format
if ~iscell(dfnl) % multi-index input
    dfnl = set_fnl_tensor(dfnl);
end

end