function [data]      = set_Fext(obj)
% Creates the data struct from the input second order force
% Structs for storing the coefficients
F_n_k = repmat(struct('coeffs',[],'ind',[]),numel(obj.fext.data(1).f_n_k),1);
data  = repmat(struct('kappa',[],'F_n_k',[]),numel(obj.fext.data),1);

% Fill the structs
for i = 1:numel(obj.fext.data)
    for j = 1:numel(obj.fext.data(i).f_n_k)
        F_n_k(j).coeffs = [obj.fext.data(i).f_n_k(j).coeffs;...
            sparse(obj.n, size(obj.fext.data(i).f_n_k(j).coeffs,2)) ];
        F_n_k(j).ind = [obj.fext.data(i).f_n_k(j).ind.';...
            sparse(obj.n, size(obj.fext.data(i).f_n_k(j).ind,1)) ].';
    end
    data(i).F_n_k = F_n_k;
    data(i).kappa = obj.fext.data(i).kappa;
end
end