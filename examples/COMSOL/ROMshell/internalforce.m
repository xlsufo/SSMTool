function L = internalforce(model, input, Null, Nullf, ud ,K0, C0)
n = length(input)/2;
u = input(1:n);
v = input(n+1:end);
u0 = Null*u + ud;
v0 = Null*v;

model.sol('sol2').setU(u0);
model.sol('sol2').setUDot(v0);
model.sol('sol2').setPNames('t');
model.sol('sol2').setPVals(1);
model.sol('sol2').createSolution;
ML = mphmatrix(model, 'sol2', 'Out', {'L','K'}, 'initmethod','sol', 'initsol', 'sol2');
Ln = Nullf'*(ML.L-ML.K*ud);

L = -Ln - K0*u - C0*v;

end