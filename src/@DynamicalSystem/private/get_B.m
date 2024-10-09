function B     = get_B(obj,obj_B)
if obj.order ==1

    if isempty(obj_B) && isempty(obj.K)
        B = speye(obj.N,obj.N);
    elseif ~isempty(obj.K) % Underlying system is second order
        B = [obj.C,    obj.M;
            obj.M,  sparse(obj.n,obj.n)];
    else
        B = obj_B;
    end

elseif obj.order == 2

    B = [obj.C,    obj.M;
        obj.M,  sparse(obj.n,obj.n)];
end
end