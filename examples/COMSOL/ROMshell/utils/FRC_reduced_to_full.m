function FRC = FRC_reduced_to_full(nt,Nonauto,contrNonauto,FRC,FRCdata,W_0,W_1,outdof)
% FRC_REDUCED_TO_FULL This function maps the forced response curve of
% equilibria/periodic orbits/two-dimensional invariant tori to the forced
% response curve of periodic orbits/two/three-dimensional invariant tori of
% the full system

% pre-processing
iNonauto = Nonauto.iNonauto;
kNonauto = Nonauto.kNonauto;
rNonauto = Nonauto.rNonauto;
mFreqs   = FRCdata.mFreqs;
m        = numel(mFreqs);


timeFRCPhysicsDomain = tic;

% check toolbox
Zout_frc  = [];
Znorm_frc = [];
Aout_frc  = [];
ZoutNorm_frc = [];
Zic_frc = [];


% Loop around a resonant mode
om   = FRC.om;
epsf = FRC.ep;

for j = 1:numel(om)
    j/numel(om)
    %% ep toolbox
    state = FRC.z(j,:);
    if contrNonauto
        W_1j = W_1{j};
    else
        W_1j = W_1;
    end
    [Aout, Zout, z_norm, Zic, ZoutNorm] = compute_full_response_2mD_ReIm(W_0, W_1j, state, epsf(j), nt, mFreqs, outdof);

    % collect output in array
    Aout_frc = [Aout_frc; Aout];
    Zout_frc = [Zout_frc; Zout];
    Znorm_frc = [Znorm_frc; z_norm];
    ZoutNorm_frc = [ZoutNorm_frc; ZoutNorm];
    Zic_frc = [Zic_frc; Zic];

end
%% 
% Record output
FRC.Aout_frc  = Aout_frc;
FRC.Zout_frc  = Zout_frc;
FRC.Znorm_frc = Znorm_frc;
FRC.ZoutNorm_frc = ZoutNorm_frc;
FRC.Zic_frc = Zic_frc;

FRC.timeFRCPhysicsDomain = toc(timeFRCPhysicsDomain);
FRC.SSMorder   = FRCdata.order;
FRC.SSMispolar = FRCdata.ispolar;

end