function F = set_F(obj,F)
% sets nonlinearity in second order form in tensor format
if ~isempty(F) && ~iscell(F) % multi-index input
    F = set_F_tensor(F);
end

end