function y = Nhan(obj,u,v)
% NHAN This function returns internal force N(u,v) in equation of motion as
% follows M\ddot{u}+N(u,\dot{u})=F(t,p). This function will be used in
% forward toolbox for finding periodic orbits
%
% Y = NHAN(OBJ,U,V)
%
% obj: coco object with dynamical system included
% u:   displacement
% v:   velocity
% 
% See also: DNDU, DNDV

y = obj.system.C*v+obj.system.K*u+obj.system.compute_fnl(u,v);

end