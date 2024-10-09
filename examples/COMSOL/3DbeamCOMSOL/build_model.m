function [model, M, C, K, f_0, Null, Nullf, ud, outdof, u0, disp_idx] = build_model(f0)
import com.comsol.model.*
import com.comsol.model.util.*

mphopen -clear;

filename_comsol = 'beam.mph';
model = mphload(filename_comsol);

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);

model.study('std1').feature('stat').set('geometricNonlinearity', true);
model.study('std2').feature('time').set('geometricNonlinearity', true);


%% get MCK from time-dependent problem
MA = mphmatrix(model ,'sol2', ...
 'Out', {'K','D','E','Null','Nullf','ud','uscale'},'initmethod','sol','initsol','zero');
M_0 = MA.E; C_0 = MA.D; K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; 

M = Nullf'*M_0*Null;
C = Nullf'*C_0*Null;
K = Nullf'*K_0*Null;

%% set boundary load to get f_ext
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 f0]);

MA = mphmatrix(model ,'sol1', ...
 'Out', {'K','Null','Nullf','ud','uscale','L'},'initmethod','sol','initsol','zero');
Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud;
f_0 = Nullf'*(MA.L-MA.K*MA.ud);
outdof = find(f_0);
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);
%% get the index for displacement
info_mesh = mphxmeshinfo(model);
for i = 1:length(info_mesh.dofs.dofnames)
    if strcmp(info_mesh.dofs.dofnames{i}, 'comp1.u')
        break;
    end
end

disp_idx = find(info_mesh.dofs.nameinds>=i-1);
u0 = mphgetu(model,'soltag','sol1');

end