function BinvA = get_BinvA(obj,obj_BinvA)
if isempty(obj_BinvA)
BinvA = [sparse(obj.n,obj.n), speye(obj.n,obj.n)
    -obj.M\obj.K,   -obj.M\obj.C];
else
    BinvA = obj_BinvA;
end
end