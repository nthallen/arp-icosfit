function CavityLength = load_cavity_length;
% CavityLength = load_cavity_lenght;
% Looks for CavityLength.mat in the current or parent directory.
dirs = { '.', '..' };
for i = 1:length(dirs)
  path = [ dirs{i} '/CavityLength.mat' ];
  if exist(path,'file')
    load(path);
    return;
  end
end
error('Could not locate CavityLength.mat');
