function input_dim = get_F_semi_input_dim(F_semi,n)
% GET_F_SEMI_INPUT_DIM
% Input dimension of semi-intrusive second order system force

try
    
    v = zeros(n,1);
    for j = 1:numel(F_semi)
        fnl_j = F_semi{j};

        input = cell(j+1,1);
        for i = 1:(j+1)
            input{i} = v;
        end
        if ~isempty(fnl_j)
            fnl_j(input)
        end
    end

    input_dim = n;

catch

    try
        v = zeros(N,1);

        for j = 1:numel(F_semi)
            fnl_j = F_semi{j};

            input = cell(j+1,1);
            for i = 1:(j+1)
                input{i} = v;
            end
            if ~isempty(fnl_j)
                fnl_j(input)
            end

        end

        input_dim = N;
    catch

        error('Please check input dimension for the nonlinearity function handle')
    end
end
