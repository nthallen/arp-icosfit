function [drift,scannum] = laser_drift(base, ptefile);
% [drift,scannum] = laser_drift(base, ptefile);
if nargin < 2 || length(ptefile) == 0
  ptefile = 'PTE.txt';
end
P = load(ptefile);
ca = P(:,1);
fn = P(:,5)*0.0198;
oldpwd = pwd;
cd(base);
ICOSconfig;
fd = load('ICOSsum.dat');
cd(oldpwd);
cb = fd(:,6);
nu_F0 = fd(:,n_base_params+n_input_params+1);
fnb = interp1(ca,fn,cb);
drift = fnb - nu_F0;
drift = drift-drift(1);
if nargout > 1
  scannum = cb;
end
