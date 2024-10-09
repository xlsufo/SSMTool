function A     = get_A(obj,obj_A)
            
if obj.order ==1

    if ~isempty(obj.K) %Underlying system is second order
        A = [-obj.K,         sparse(obj.n,obj.n);
            sparse(obj.n,obj.n),   obj.M];
    else
        A = obj_A;
    end
elseif obj.order == 2
    A = [-obj.K,         sparse(obj.n,obj.n);
        sparse(obj.n,obj.n),   obj.M];
end

% For 
end