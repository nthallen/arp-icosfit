function fsr = load_FSR;
% fsr = load_FSR;
% Looks for EtalonFSR.mat in the current or parent directory.
dirs = { '.', '..' };
for i = 1:length(dirs)
  path = [ dirs{i} '/EtalonFSR.mat' ];
  if exist(path,'file')
    load(path);
    return;
  end
end
fsr = 0.0198;
