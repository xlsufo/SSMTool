function add_forcing(obj,f,varargin)
if ~isfield(f, 'data') % new format
    switch obj.order
        case 1
            nn = size(f,1);
            Kappas = varargin{1};
            data(1).kappa = Kappas(1);
            data(2).kappa = Kappas(2);
            if nn == obj.N
                data(1).F_n_k(1).coeffs = f(:,1);
                data(1).F_n_k(1).ind    = sparse(1,obj.N);
                data(2).F_n_k(1).coeffs = f(:,2);
                data(2).F_n_k(1).ind    = sparse(1,obj.N);
            elseif nn == obj.N/2 % second order system, first order computation
                data(1).F_n_k(1).coeffs = [f(:,1);sparse(nn,1)];
                data(1).F_n_k(1).ind    = sparse(1,obj.N);
                data(2).F_n_k(1).coeffs = [f(:,2);sparse(nn,1)];
                data(2).F_n_k(1).ind    = sparse(1,obj.N);
            else
                error('Forcing vector has wrong dimension.')
            end
            f_ext.data = data;

        case 2
            nn = size(f,1);
            Kappas = varargin{1};
            data(1).kappa = Kappas(1);
            data(2).kappa = Kappas(2);
            data(1).f_n_k(1).coeffs = f(:,1);
            data(1).f_n_k(1).ind    = sparse(1,nn);
            data(2).f_n_k(1).coeffs = f(:,2);
            data(2).f_n_k(1).ind    = sparse(1,nn);
            f_ext.data = data;
    end

    if numel(varargin)>1
        f_ext.epsilon = varargin{2};
    else
        f_ext.epsilon = 1;
    end
else
    f_ext = f;
end

switch obj.order

    case 1

        obj.Fext = f_ext;

        if isfield(f_ext,'epsilon')
            obj.Fext.epsilon = f_ext.epsilon;

        elseif nargin == 3
            obj.Fext.epsilon = varargin{1};

        else
            obj.Fext.epsilon = 1;

        end

    case 2
        obj.fext.data = f_ext.data;

        if isfield(f_ext,'epsilon')
            obj.fext.epsilon = f_ext.epsilon;

        elseif nargin == 3
            obj.fext.epsilon = varargin{1};

        else
            obj.fext.epsilon = 1;

        end
end
end