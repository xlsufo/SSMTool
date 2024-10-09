function prob = cocoSet(obj, prob)
%COCOSET This function present settings for coco

opts    = obj.Options;
default = cocoOptions();
% settings for continuation
if opts.NPR~=default.NPR
    prob = coco_set(prob, 'cont', 'NPR', opts.NPR);
end
if opts.NSV~=default.NSV
    prob = coco_set(prob, 'cont', 'NSV', opts.NSV);
end
if opts.NAdapt~=default.NAdapt
    prob = coco_set(prob, 'cont', 'NAdapt', opts.NAdapt);
end
if opts.h0~=default.h0
    if strcmp(obj.atlasAlg,'1d')
        prob = coco_set(prob, 'cont', 'h0', opts.h0);
    else
        prob = coco_set(prob, 'cont', 'atlas', 'kd', 'R', opts.h0);
    end
end
if opts.h_max~=default.h_max
    if strcmp(obj.atlasAlg,'1d')
        prob = coco_set(prob, 'cont', 'h_max', opts.h_max);
    else
        prob = coco_set(prob, 'cont', 'atlas', 'kd', 'R_max', opts.h_max);
    end
end
if opts.h_min~=default.h_min
    if strcmp(obj.atlasAlg,'1d')
        prob = coco_set(prob, 'cont', 'h_min', opts.h_min);
    else
        prob = coco_set(prob, 'cont', 'atlas', 'kd', 'R_min', opts.h_min); 
    end
end
if opts.h_fac_max~=default.h_fac_max
    if strcmp(obj.atlasAlg,'1d')
        prob = coco_set(prob, 'cont', 'h_fac_max', opts.h_fac_max);
    else
        prob = coco_set(prob, 'cont', 'atlas', 'kd', 'R_fac_max', opts.h_fac_max); 
    end
end
if opts.h_fac_min~=default.h_fac_min
    if strcmp(obj.atlasAlg,'1d')
        prob = coco_set(prob, 'cont', 'h_fac_min', opts.h_fac_min);
    else
        prob = coco_set(prob, 'cont', 'atlas', 'kd', 'R_fac_min', opts.h_fac_min);
    end
end
if opts.al_max~=default.al_max
    if strcmp(obj.atlasAlg,'1d')
        prob = coco_set(prob, 'cont', 'al_max', opts.al_max);
    else
        prob = coco_set(prob, 'cont', 'atlas', 'kd', 'almax', opts.al_max);
    end
end
if opts.MaxRes~=default.MaxRes
    prob = coco_set(prob, 'cont', 'MaxRes', opts.MaxRes);
end
if opts.bi_direct~=default.bi_direct
    prob = coco_set(prob, 'cont', 'bi_direct', opts.bi_direct);
end
if opts.PtMX~=default.PtMX
    prob = coco_set(prob, 'cont', 'PtMX', opts.PtMX);
end
if opts.theta~=default.theta && strcmp(obj.atlasAlg,'kd')
    prob = coco_set(prob, 'cont', 'atlas', 'kd', 'theta', opts.theta);
end
% settings for correction
if opts.ItMX~=default.ItMX
    prob = coco_set(prob, 'corr', 'ItMX', opts.ItMX);
end
if opts.TOL~=default.TOL
    prob = coco_set(prob, 'corr', 'TOL', opts.TOL);
end
% settings for collocation
if opts.NTST~=default.NTST
    prob = coco_set(prob, 'coll', 'NTST', opts.NTST);
end
if opts.NCOL~=default.NCOL
    prob = coco_set(prob, 'coll', 'NCOL', opts.NCOL);
end
if opts.MXCL~=default.MXCL
    prob = coco_set(prob, 'coll', 'MXCL', opts.MXCL);
end
% settings for forward
if opts.IntItMX~=default.IntItMX
    ode_opts.ItMX = opts.IntItMX;
    prob = coco_set(prob, 'forward', 'ode_opts', ode_opts);
end
if opts.RelTol~=default.RelTol
    ode_opts.RelTol = opts.RelTol;
    prob = coco_set(prob, 'forward', 'ode_opts', ode_opts);
end
if opts.Nsteps~=default.Nsteps
    ode_opts.Nsteps = opts.Nsteps;
    prob = coco_set(prob, 'forward', 'ode_opts', ode_opts);
end
if opts.alpha~=default.alpha
    ode_opts.alpha = opts.alpha;
    prob = coco_set(prob, 'forward', 'ode_opts', ode_opts);
end
if opts.rhoinf~=default.rhoinf
    ode_opts.rhoinf = opts.rhoinf;
    prob = coco_set(prob, 'forward', 'ode_opts', ode_opts);
end
if opts.neigs~=default.neigs
    prob = coco_set(prob, 'po', 'neigs', opts.neigs);
end
end