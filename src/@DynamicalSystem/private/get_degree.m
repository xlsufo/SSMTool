function degree  = get_degree(obj)
% GET_DEGREE
% Get method

degree = 0;
if ~isempty(obj.A)
    degree = length(obj.F);
end

end