function WaveSpecs = load_waves(errok);
% WaveSpecs = load_waves;
% Find waves.m, run it and return the structure it defines.
if nargin < 1
  errok = 0;
end
run = getrun(1);
cr_cfg = load_cr_cfg;
dirs = { '.', [ cr_cfg.Matlab_CD_Path cr_cfg.HomeDir '/anal/' run ], ...
    [ cr_cfg.Matlab_CD_Path cr_cfg.HomeDir filesep run ], ...
    [ cr_cfg.Matlab_CD_Path cr_cfg.HomeDir filesep run '/Base' ] };
for j = 1:length(dirs);
  path = [ dirs{j} '/waves.m' ];
  if exist(path,'file')
    % I did this before with 'run()', but run doesn't deal
    % with forward slashes.
    savedir = pwd;
    cd(dirs{j});
    waves;
    cd(savedir);
    return;
  end
end
if errok
  WaveSpecs = [];
else
  error([ 'Unable to locate waves.m' ]);
end
