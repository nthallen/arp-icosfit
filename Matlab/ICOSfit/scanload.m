function fe = scanload( ScanNum, base, binary );
% f = scanload( CPCI14 [, base[, binary]] )
% binary defaults to TRUE
% base defaults to Config File ScanDir
% Looks for the file either directly under base or
% prefixes base with 'E:/home/CR/<run>/CPCI/' to look
% on a CD (assuming the CD is at E:)
ICOSfit_cfg = load_ICOSfit_cfg;
if nargin < 3
  binary = 1;
  if nargin < 2
    base = ICOSfit_cfg.ScanDir;
  end
end
ifile = mlf_path(base,ScanNum,'.dat');
run = getrun(1);
iifile = ifile;
if ~exist( iifile, 'file')
  iifile = [ ICOSfit_cfg.Matlab_Path filesep run filesep ifile ];
end
if binary == 0
  fe = load(iifile);
else
  fid = fopen(iifile, 'r');
  if fid < 0
    error([ 'Unable to locate file "' ifile '"' ]);
  end
  [ dim, count ] = fread(fid,2,'uint32');
  fe = fread(fid,dim','float32');
  fclose(fid);
end
