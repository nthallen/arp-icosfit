function PTE_gen( cpci, PTEfile, PTparams );
% PTE_gen( [ cpci [, PTEfile [, PTparams ] ] ] );
%  Make up numbers for a tuning rate in the absence of an etalon
%  If we don't have P and T, make them up too.
%  For cpci, use cpci_list
% PTparams correspond to columns 4:end of PTE.txt
if nargin < 3
  PTparams = ...
    [200 17.91724 34.33037 1.489716 -9.231663 0.09478366 -14.00556 0.2113355 ];
elseif any(size(PTparams) ~= [1 8])
  error('PTparams input must be 1x8');
end
if nargin < 2
  PTEfile = 'PTE.txt';
end
if nargin < 1
  cpci = cpci_list;
end

if exist('PT.mat','file')
  PT = load('PT.mat');
else
  PT = [];
end
if isfield(PT,'CPCI14')
  if min(cpci) < min(PT.CPCI14) | max(cpci) > max(PT.CPCI14)
    warning('CPCI14 files exceed PT file indices');
  end
  v = find(diff(PT.CPCI14)>0)+1;
end
if isfield(PT,'CPCI14') && isfield(PT,'CellP')
  P = interp1(PT.CPCI14(v),PT.CellP(v),cpci,'nearest');
else
  warning('Assuming 40 Torr');
  P = 40*ones(size(cpci));
end
if isfield(PT,'CPCI14') && isfield(PT,'Tavg')
  T = interp1(PT.CPCI14(v),PT.Tavg(v),cpci,'nearest');
else
  warning('Assuming 25 Celcius');
  T = (25+273.15)*ones(size(cpci));
end
ofp = fopen(PTEfile, 'a');
%fprintf( ofp, '%d %.2f %.1f %d %.7g %.7g %.7g %.7g %.7g %.7g %.7g\n', ...
%  0, P(1), T(1), PTparams );
for i=1:length(cpci)
  fprintf( ofp, '%d %.2f %.1f %d %.7g %.7g %.7g %.7g %.7g %.7g %.7g\n', ...
    cpci(i), P(i), T(i), PTparams );
end
fclose(ofp);

      