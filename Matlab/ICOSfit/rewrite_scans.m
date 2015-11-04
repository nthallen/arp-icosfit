function rewrite_scans(ibase,obase,M)
%This program is useful if spectra are inverted. It is better to fix it in hardware. 
% Should be customized for each axis and probably a specific version should be made 
% and copied into local directories.
if nargin < 1 || isempty(ibase)
   ibase = 'SSP';
end
if nargin < 2 || isempty(obase)
   obase = 'SSPo';
end
if nargin < 3 || isempty(M)
    M = [ 0 0 1; 0 1 0; 1 0 0];
end
    
index = 1;
while 1
  pi = mlf_path(ibase,index);
  [fi,hdr] = loadbin(pi);
  if isempty(fi); break; end
  po = mlf_path(obase,index);
  fo = fi * M;
  mlf_mkdir(obase,index);
  writebin( po, fo, hdr );
  index = index+1;
end
