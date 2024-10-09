function y = Fext(obj,t,p)
% FEXT This function returns external force F in equation of motion as
% follows M\ddot{u}+N(u,\dot{u})=F(t,p). This function will be used in
% forward toolbox for finding periodic orbits
%
% Y = FEXT(OBJ,T,P)
%
% obj: coco object with dynamical system included
% t:   time
% p:   problem parameter (omega,epsilon)
% 
% See also: DFEXTDP

om = p(1);
ep = p(2);
fext = obj.system.fext.data.f_n_k; kapa = obj.system.fext.data.kappa;
fext_coeffs = fext.coeffs;
fext_harm   = ep*cos(kapa(1)*om*t);
y = 2*fext_coeffs*fext_harm;

end

