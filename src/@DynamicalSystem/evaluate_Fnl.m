function Fnl = evaluate_Fnl(obj,z)
% EVALUATE_FNL this function evaluates the nonlinearity at input state z
% in a dynamical system object in the first-order form

if strcmp(obj.Options.Intrusion,'none')
    switch obj.order
        case 1
            Fnl_handle = obj.F_non;
        case 2
            if obj.nl_input_dim == obj.N
                Fnl_handle = @(v) [-obj.fnl_non(v); sparse(obj.n,1)];
            else
                Fnl_handle = @(v) [-obj.fnl_non(v(1:obj.n)); sparse(obj.n,1)];
            end
    end
    Fnl = zeros(size(z));
    for i= 1:size(z,2)
        % If input function is vectoriseable, then this could be vectorised
        Fnl(:,i) = Fnl_handle(z(:,i));

    end

else
    switch obj.order
        case 1
            Fnl = zeros(obj.N,1);
            degree = length(obj.F);
            if iscell(obj.F)
                % case tensor nonlinearity
                for j = 2:degree
                    % expand_tensor does not support vectorization
                    % Fnl = Fnl + expand_tensor(obj.F{j},z);
                    Fj  = tensor_to_multi_index(obj.F{j});
                    Fnl = Fnl + expand_multiindex(Fj,z);
                end
            else   % case multi-index nonlinearity
                for j = 2:degree
                    Fnl = Fnl + expand_multiindex(obj.F(j),z);
                end
            end

        case 2
            [x, xd] = deal(z(1:obj.n,:),z(obj.n+1:obj.N,:));
            nt = size(x,2);
            fnl = obj.compute_fnl(x,xd);
            Fnl = [-fnl;
                sparse(obj.n,nt)];
    end
end