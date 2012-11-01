function base = find_scans_dir( base_in, ICOSfit );
% base = find_scan_dir(base_in [, ICOSfit]);
% Finds the appropriate base directory for log files
% If the ICOSfit argument is provided and is non-zero, the
% path appropriate for use by ICOSfit is returned.
% (see load_ICOSfit_cfg for details).
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
    ICOSfit_cfg = load_ICOSfit_cfg;
    base = ICOSfit_cfg.ScanDir;
    if ~exist ( base, 'dir' )
    base = [ ICOSfit_cfg.Matlab_Path '/' getrun(1) ...
       '/' ICOSfit_cfg.ScanDir ];
    end
    if ~exist( base, 'dir' )
      error( [ 'Unable to locate appropriate base for ' getrun(0) ] );
    end
    if nargin > 1 && ICOSfit
      base = [ ICOSfit_cfg.ICOSfit_Path '/' getrun(1) ...
         '/' ICOSfit_cfg.ScanDir ];
    end
end
