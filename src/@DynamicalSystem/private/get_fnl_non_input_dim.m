function input_dim = get_fnl_non_input_dim(fnl_non,n,N)
% GET_FNL_NON_INPUT_DIM
% Input dimension of intrusive second order force

try
    
    v = zeros(n,1);
    fnl_non(v);
    input_dim = n;

catch

    try
        v = zeros(N,1);
        fnl_non(v);
        input_dim = N;
    catch

        error('Please check input dimension for the nonlinearity function handle')
    end
end


end