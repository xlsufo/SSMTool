function y = dNdu(obj,u,v)
% DNDU This function returns the derivative of internal force N(u,v) with
% respect to its first argument. EOM is given by M\ddot{u}+N(u,\dot{u})=F(t,p). 
% This function will be used in forward toolbox for finding periodic orbits
%
% Y = DNDU(OBJ,U,V)
%
% obj: coco object with dynamical system included
% u:   displacement
% v:   velocity
% 
% See also: NHAN, DNDV

y = obj.system.K+obj.system.compute_dfnldx(u,v);

% Jx = coco_ezDFDX('f(x,p)', @obj.Nhan, u, v);

end