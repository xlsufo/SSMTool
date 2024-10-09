function nl_input_dim = get_nl_input_dim(obj,obj_nl_input_dim)
% GET_NL_INPUT_DIM
% Get method

if ~isempty(obj_nl_input_dim)
    nl_input_dim = obj_nl_input_dim;

else
    if ~isempty(obj.fnl)
        nl_input_dim = get_fnl_input_dim(obj.fnl);

    elseif ~isempty(obj.F)
        nl_input_dim = get_F_input_dim(obj.F);

    elseif ~isempty(obj.fnl_non)
        nl_input_dim = get_fnl_non_input_dim(obj.fnl_non,obj.n,obj.N);

    elseif ~isempty(obj.F_non)
        nl_input_dim = get_F_non_input_dim(obj.F_non,obj.n,obj.N);

    elseif ~isempty(obj.fnl_semi)
        nl_input_dim = get_fnl_semi_input_dim(obj.fnl_semi,obj.n,obj.N);

    elseif ~isempty(obj.F_semi)
        nl_input_dim = get_F_semi_input_dim(obj.F_semi,obj.n,obj.N);
    else
        error("No input dimension for intrusive nonlinearity set - please check F / fnl")
    end

    obj.nl_input_dim = nl_input_dim;
end

end