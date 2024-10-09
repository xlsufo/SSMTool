function [F_non]     = set_F_non(obj)

F_non = @(v) [-obj.fnl_non(v); sparse(obj.n,1)];

end