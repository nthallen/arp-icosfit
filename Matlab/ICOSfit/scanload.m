function fe = scanload( ScanNum, base, binary )
% f = scanload( ScanNum [, base[, binary]] )
% binary defaults to TRUE
% base defaults to Config File ScanDir
% Looks for the file either directly under base or
% prefixes base with 'E:/home/CR/<run>/ScanDir/' to look
% on a CD (assuming the CD is at E:)
ICOSfit_cfg = load_ICOSfit_cfg;
if nargin < 3
  binary = 1;
  if nargin < 2
    base = '';
  end
end
base = find_scans_dir(base);
ifile = mlf_path(base,ScanNum,'.dat');
if binary == 0
  fe = load(ifile);
else
  fe = loadbin(ifile);
end
