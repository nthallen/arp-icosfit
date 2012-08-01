function fe = scanload( CPCI14, base, binary );
% f = scanload( CPCI14 [, base[, binary]] )
% binary defaults to TRUE
% base defaults to CPCI14
% Looks for the file either directly under base or
% prefixes base with 'E:/home/CR/<run>/CPCI/' to look
% on a CD (assuming the CD is at E:)
cr_cfg = load_cr_cfg;
if nargin < 3
  binary = 1;
  if nargin < 2
    base = 'CPCI14';
  end
end
ifile = mlf_path(base,CPCI14,'.dat');
run = getrun(1);
iifile = ifile;
if ~exist( iifile, 'file')
  iifile = [ cr_cfg.Matlab_CD_Path cr_cfg.HomeDir filesep run cr_cfg.CPCI14link ifile ];
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
