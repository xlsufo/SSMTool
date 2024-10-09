function input_dim = get_F_non_input_dim(F_non,n,N)
% GET_F_NON_INPUT_DIM
% Input dimension of non-intrusive first order force


try
    
    v = zeros(n,1);
    F_non(v);
    input_dim = n;

catch

    try
        v = zeros(N,1);
        F_non(v);
        input_dim = N;
    catch

        error('Please check input dimension for the nonlinearity function handle')
    end
end


end