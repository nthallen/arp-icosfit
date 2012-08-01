function base = find_scan_dir( base_in, ICOSfit );
% base = find_scan_dir(base_in [, ICOSfit]);
% Finds the appropriate base directory for CPCI14 log files
% If the ICOSfit argument is provided and is non-zero, the
% path appropriate for use by ICOSfit is returned.
% (see load_cr_cfg for details).
if nargin < 1
  base_in = '';
end
if length(base_in)
  if exist(base_in,'dir')
    base = base_in;
  else
    error(sprintf('Input base "%s" is not a directory', base_in));
  end
else
  if exist('CPCI14','dir')
    base = 'CPCI14';
  else
    cr_cfg = load_cr_cfg;
    base = [ cr_cfg.Matlab_CD_Path cr_cfg.HomeDir '/' getrun(1) ...
      cr_cfg.CPCI14link 'CPCI14' ];
    if ~exist( base, 'dir' )
      error( [ 'Unable to locate appropriate base for ' getrun(0) ] );
    end
    if nargin > 1 && ICOSfit
      base = [ cr_cfg.ICOSfit_CD_Path cr_cfg.HomeDir '/' getrun(1) ...
        cr_cfg.CPCI14link 'CPCI14' ];
    end
  end
end
