function intrusive_input_dim = get_fnl_input_dim(fnl)
% GET_FNL_INPUT_DIM
%
d = length(fnl) ;

j = 1;
intrusive_input_dim = [];
while isempty(intrusive_input_dim) || j <= d
    
    if ~isempty(fnl{j}) && nnz(fnl{j})
        
    intrusive_input_dim = size(fnl{j},2);
        
    end

    j=j+1;
    
end

if isempty(intrusive_input_dim)
    error("Failed to set input dimension of nonlinearity")
end
end
