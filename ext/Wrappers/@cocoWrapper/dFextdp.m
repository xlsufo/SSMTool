function y = dFextdp(obj,t,p)
% FEXT This function returns the derivative of external force F with
% respect to problem parameters p. EOM: M\ddot{u}+N(u,\dot{u})=F(t,p). 
% This function will be used in forward toolbox for finding periodic orbits
%
% Y = DFEXTDP(OBJ,T,P)
%
% obj: coco object with dynamical system included
% t:   time
% p:   problem parameter (omega,epsilon)
% 
% See also: FEXT

om = p(1);
ep = p(2);
fext = obj.system.fext.data.f_n_k; 
kapa = obj.system.fext.data.kappa;
fext_coeffs = fext.coeffs;
fext_harm_ep = cos(kapa(1)*om*t);
fext_harm_om = -kapa(1)*t*ep*sin(kapa(1)*om*t);
y = 2*fext_coeffs*[fext_harm_om fext_harm_ep];

% Jp = coco_ezDFDP('f(x,p)', @obj.Fext, t, p);

end