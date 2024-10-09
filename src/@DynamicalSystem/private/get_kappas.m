function kappas  = get_kappas(obj)
% GET_KAPPAS
% Get method

%kappas stored in rows
sz_kappa = size(obj.Fext.data(1).kappa,2);
kappas = reshape([obj.Fext.data.kappa],sz_kappa,[]).';

end