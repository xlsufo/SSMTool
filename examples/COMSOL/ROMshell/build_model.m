function [model, M, C, K, Null, Nullf, ud, outdof] = build_model()
import com.comsol.model.*
import com.comsol.model.util.*

mphopen -clear;

filename_comsol = 'plate.mph';
model = mphload(filename_comsol);

%%% get MCK from time-dependent problem
MA = mphmatrix(model ,'sol2', ...
 'Out', {'K','D','E','Null','Nullf','ud','uscale'},'initmethod','sol','initsol','zero');
M_0 = MA.E; C_0 = MA.D; K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud;

M = Nullf'*M_0*Null;
C = Nullf'*C_0*Null;
K = Nullf'*K_0*Null;
%% get the outdof in the middle of the plate
info_mesh = mphxmeshinfo(model);
coords = info_mesh.dofs.coords;
[~,midpoint_idx] = min(vecnorm(coords-[0;0;0.02913]));  % find the index correspond to mid point with given coordinates
midpoint_idx = midpoint_idx + 5;

n = size(Null,1);
Out = zeros(n,1); Out(midpoint_idx) = 1; out = Null'*Out; outdof = find(out>0);
end