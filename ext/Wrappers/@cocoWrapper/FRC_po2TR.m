function bd = FRC_po2TR(obj, run, lab, parRange, varargin)
% FRC_PO2TR: This function performs the continuation of Neimark-Sacker (TR)
% bifurcation periodic orbits using po-toolbox of coco. In such a
% continuation, both forcing frequency omega and its amplitude are free to
% change. 
% 
% BD = FRC_PO2TR(OBJ,RUN,LAB,PARRANGE,VARARGIN)
%
% run:      runid of continuation for saved solution
% lab:      label of continuation for saved solution, which must be the
%           label of a saddle-node point
% parRange: continuation domain of parameters. It is of the form
%           {[om1,om2],[f1,f2]}, where [om1,om2] and [f1,f2] specify the
%           continuation domain of excitation frequency and amplitude
%           respectively. You can give empty array and then no domain is
%           specified, e.g., {[],[f1,f2]} only presents the domain of
%           forcing amplitude
% varargin: ['saveICs'] flag for saving the initial point in trajectory
%
% See also: EXTRACT_FRC, FRC_TR2TOR

prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = ode_po2TR(prob, '', coco_get_id(run, 'FRC'), lab);

[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
omData = struct();
omData.periodsRatio = obj.periodsRatio;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);

% track amplitude of outdof
outdof = obj.outdof;
ampdata.dof  = outdof;
ampdata.zdim = obj.system.N;
numoutdof = numel(outdof);
ampNames = cell(1, numel(numoutdof));
for k = 1:numel(outdof)
   ampNames{k} = strcat('amp',num2str(outdof(k))); 
end
prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
    'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);

cont_args = {1, [{'omega'} {'eps'} {'po.period'} ampNames(:)'], parRange};

runid = coco_get_id(obj.Options.dir_name, 'TR');
fprintf('\n Run=''%s'': Continue primary family of TR bifurcation periodic orbits.\n', ...
  runid);
bd    = coco(prob, runid, [], cont_args{:});

end
