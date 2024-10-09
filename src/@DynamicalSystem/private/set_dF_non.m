function [dF_non]    = set_dF_non(obj)

dF_non = @(v) [-obj.dfnl_non(v), sparse(obj.n,obj.n); sparse(obj.n,2*obj.n)];

end
