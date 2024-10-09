function bd = FRC_isol2tor(obj, omegaRange, t0, x0, p0)
% FRC_ISOL2TR: This function performs the continuation of Neimark-Sacker (TR)
% bifurcation periodic orbits using po-toolbox of coco. In such a
% continuation, both forcing frequency omega and its amplitude are free to
% change. 
% 
% BD = FRC_ISOL2TOR(OBJ,RUN,LAB,OMEGARANGE,NSEG,ORIENT,EPSILON)
%
% oid:        runid of current continuation
% run:        runid of continuation for saved solution
% lab:        label of continuation for saved solution, which must be the
%             label of a saddle-node point
% omegaRange: continuation domain of forcing frequency
% t0:         time points of a mesh
% x0:         a collection of trajectories in the format (nt, dim, nseg)
% p0:         problem parameters
%
% See also: EXTRACT_FRC, FRC_TR2TOR

prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'tor', 'autonomous', false);
prob = coco_set(prob, 'tor', 'Om2idx', 1);

obj.fnlTensor2Multi();
odedata.fnl = obj.multiFnl;
odefun = @(t,x,p) obj.ode_het(t,x,p,odedata);
funcs = {odefun};

tor_args = [funcs(:)', {t0}, {x0}, {{'omega','eps','om1','om2','varrho'}}, {p0}];
prob = ode_isol2tor(prob, '',tor_args{:});

cont_args = {1, {'varrho' 'omega' 'om1' 'om2' 'eps'}, {[],omegaRange}};
runid = coco_get_id(obj.Options.dir_name, 'isol2tor');
fprintf('\n Run=''%s'': Continue tori born from a TR periodic orbit.\n', ...
  runid);
bd    = coco(prob, runid, [], cont_args{:});

end
