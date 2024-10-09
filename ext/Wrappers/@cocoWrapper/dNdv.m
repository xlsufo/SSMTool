function y = dNdv(obj,u,v)
% DNDV This function returns the derivative of internal force N(u,v) with
% respect to its second argument. EOM is given by M\ddot{u}+N(u,\dot{u})=F(t,p). 
% This function will be used in forward toolbox for finding periodic orbits
%
% Y = DNDV(OBJ,U,V)
%
% obj: coco object with dynamical system included
% u:   displacement
% v:   velocity
% 
% See also: NHAN, DNDU

y = obj.system.C+obj.system.compute_dfnldxd(u,v);

end