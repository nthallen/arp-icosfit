function cr_cfg = load_cr_cfg(use_default);
% cr_cfg = load_cr_cfg;
% Checks for the existence of a CR_Config script and complains
% if one cannot be found. Below is a sample version of a
% CR_Config function.
%
%  function cr_cfg = CR_Config;
%  % CR_Config defines local configuration and should be located
%  % in a local configuration directory (as opposed to a directory
%  % of shared scripts).
%  cr_cfg.Matlab_CD_Path = 'E:';
%  cr_cfg.ICOSfit_CD_Path = '/cd';
%  cr_cfg.HomeDir = '/home/CR';
if nargin < 1
  use_default = 0;
end
dirs = { '.', '..', '../..' };
loaded = 0;
for i=1:length(dirs)
  fname = [ dirs{i} '/CR_Config.m' ];
  if exist(fname, 'file')
    cur = cd;
    cd( dirs{i} );
    cr_cfg = CR_Config;
    cd( cur );
    loaded = 1;
    break;
  end
end
if ~loaded
  if ispc
    cr_cfg.Matlab_CD_Path = 'E:';
    cr_cfg.ICOSfit_CD_Path = '/cd';
  else
    cr_cfg.Matlab_CD_Path = 'cd';
    cr_cfg.ICOSfit_CD_Path = 'cd';
  end

end
if ~isfield(cr_cfg,'HomeDir')
  cr_cfg.HomeDir = '/home/CR';
end
if ~isfield(cr_cfg,'CPCI14link')
  cr_cfg.CPCI14link = '/CPCI/';
end
if ~loaded && ~use_default
  cr_cfg = edit_cr_cfg( 'Config', cr_cfg );
end
