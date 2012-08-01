function [ base, etln, bkgd ] = loadscans(basepath, cpci14, x, zero);
% [ icos, etln, bkgd ] = loadscans(basepath, cpci14[, x[, zero]]);
% Extracts a number of ICOS scans from the raw files.
% base has length(x) rows and length(cpci14) columns
% x defaults to the entire scan;
% etln and bkgd are an optional outputs.
% The optional zero input if non-zero indicates that
% the mean of samples [1:wv.TzSamples] should be subtracted.
% Setting zero to 0 disables this subtraction.
if nargin < 3
  x = [];
end
if nargin < 4
  zero = 1;
elseif zero ~= 0
  zero = 1;
end
basepath = find_CPCI14_dir(basepath);

nscans = length(cpci14);
wv = waves_used(cpci14);
if length(wv) > 1
  error('More than one waveform used');
end
x0 = [5:wv.TzSamples];
if isempty(x0)
  x0 = 1;
  zero = 0;
end
if length(x)
  base = zeros(length(x), nscans);
  if nargout > 1
    etln = zeros(length(x), nscans );
  end
end
for i=1:nscans
  path=mlf_path(basepath,cpci14(i));
  if ~exist(path, 'file')
    error(['File not found: ' path]);
  end
  fe = loadbin(path);
  if length(x) == 0
    x = 1:size(fe,1);
    base = zeros(length(x),nscans);
    if nargout > 1
      etln = zeros(length(x), nscans );
    end
    if nargout > 2
        bkgd = zeros(length(x), nscans );
    end
  end
  if size(fe,1) < max(x)
    error(['Bad input at CPCI14 ' num2str(cpci14(i)) ', ' path ]);
  end
  base(:,i) = fe(x,1) - zero * mean(fe(x0,1));
  if nargout > 1
    etln(:,i) = fe(x,2) - zero * mean(fe(x0,2));
  end
  if nargout > 2
      bkgd(:,i) = fe(x,3) - zero * mean(fe(x0,3));
  end
end