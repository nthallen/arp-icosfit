function WaveSpecs = load_waves(errok)
% WaveSpecs = load_waves;
% Find waves.m, run it and return the structure it defines.
if nargin < 1
  errok = 0;
end
rundir = getrun(1);
ICOSfit_cfg = load_ICOSfit_cfg;
dirs = { '.', [ ICOSfit_cfg.Matlab_Path '/anal/' rundir ], ...
    [ ICOSfit_cfg.Matlab_Path filesep rundir ], ...
    [ ICOSfit_cfg.Matlab_Path filesep rundir '/Base' ] };
for j = 1:length(dirs);
  path = [ dirs{j} '/' ICOSfit_cfg.WavesFile ];
  if exist(path,'file')
    % I did this before with 'run()', but run doesn't deal
    % with forward slashes.
    savedir = pwd;
    cd(dirs{j});
     eval(['run ' ICOSfit_cfg.WavesFile]);
    cd(savedir);
    return;
  end
end
if errok
  WaveSpecs = [];
else
  error([ 'Unable to locate waves.m' ]);
end
