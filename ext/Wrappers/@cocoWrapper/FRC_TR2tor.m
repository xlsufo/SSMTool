function bd = FRC_TR2tor(obj, run, lab, omegaRange, nseg, epsilon,varargin)
% FRC_PO2TR: This function performs the continuation of Neimark-Sacker (TR)
% bifurcation periodic orbits using po-toolbox of coco. In such a
% continuation, both forcing frequency omega and its amplitude are free to
% change. 
% 
% BD = FRC_TR2TOR(OBJ,RUN,LAB,OMEGARANGE,NSEG,EPSILON)
%
% oid:        runid of current continuation
% run:        runid of continuation for saved solution
% lab:        label of continuation for saved solution, which must be the
%             label of a saddle-node point
% omegaRange: continuation domain of forcing frequency
% nseg:       number of segments (2nseg+1) - 10 or bigger in general
% epsilon     perturbation amount for initial torus - 1e-4 or smaller in
%             general
% varargin:   samples for omega (used in coco_add_event)
%
% See also: EXTRACT_FRC, FRC_TR2TOR

prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = coco_set(prob, 'tor', 'autonomous', false);
prob = coco_set(prob, 'tor', 'Om2idx', 1);
prob = ode_TR2tor(prob, '', run, lab, nseg, 'pos', epsilon);

if ~isempty(varargin)
    omSamp = varargin{1};
    prob   = coco_add_event(prob, 'UZ', 'omega', omSamp);
end
cont_args = {1, {'varrho' 'omega' 'om1' 'om2' 'eps'}, {[],omegaRange}};

runid = coco_get_id(obj.Options.dir_name, 'tor');
fprintf('\n Run=''%s'': Continue tori born from a TR periodic orbit.\n', ...
  runid);
bd    = coco(prob, runid, [], cont_args{:});

end
