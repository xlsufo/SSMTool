function Kn = internalforceJacobian(model,u,Null,Nullf,ud, u0, disp_idx)
u0(disp_idx) = Null*u + ud;
model.sol('sol1').setU(u0);
model.sol('sol1').setPVals(0);
model.sol('sol1').createSolution;
ML = mphmatrix(model, 'sol1', 'Out', {'K'}, 'initmethod','sol', 'initsol', 'sol1');
Kn = Nullf'*ML.K*Null;

end