function intrusive_input_dim = get_F_input_dim(F)
% GET_F_INPUT_DIM

d = length(F) ;

j = 1;

intrusive_input_dim = [];
while isempty(intrusive_input_dim) || j <= (d)
    
    if ~isempty(F{j}) && nnz(F{j})
        
    intrusive_input_dim = size(F{j},2);
        
    end

    j=j+1;
    
end

if isempty(intrusive_input_dim)
    error("Failed to set input dimension of nonlinearity")
end
end
