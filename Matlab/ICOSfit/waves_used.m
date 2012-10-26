function [ wvsused, ranges ] = waves_used(cpci14)
% waves_used( [cpci14] );
%   Lists waveforms used in the current run.
% wvs = waves_used([cpci14]);
%   Returns a struct array of waveform definitions
% [wvs,ranges] = waves_used([cpci14]);
%   Returns two struct arrays of the same length
%   ranges has 'wvno' and 'ranges' fields.
%   The ranges field is an n x 2 matrix where the two columns
%   are starting and ending CPCI numbers for each region
WaveSpecs = load_waves;
PT = load_mat_files('PT');
dcpi = find(diff(PT.CPCI14)>0)+1; % index of new cpci numbers
if nargin > 0
  dcpi = dcpi(ceil(interp1(PT.CPCI14(dcpi),1:length(dcpi),cpci14,'linear','extrap')));
  dcpi = dcpi(find(~isnan(dcpi)));
else
  cpci14 = PT.CPCI14;
end
wvs = PT.QCLI_Wave(dcpi);
wvnos = unique(wvs)';
if nargout > 0
  wvsused = WaveSpecs(wvnos+1);
end
if nargout == 1; return; end
if nargout >= 2; ranges = []; end
for wvno=wvnos
  dw = diff([-1;wvs;-1] == wvno);
  wvbeg = PT.CPCI14(dcpi((find(dw > 0)))-1)+1;
  wvend = PT.CPCI14(dcpi(find(dw < 0)-1));
  if length(wvbeg) ~= length(wvend)
    error('Different lengths');
  end
  wvrng = struct('wvno',wvno,'ranges', [ wvbeg wvend ]);
  if nargout >= 2
    if isempty(ranges)
      ranges = wvrng;
    else
      ranges(end+1) = wvrng;
    end
  else
    fprintf( 1, '%d %s:', wvno, WaveSpecs(wvno+1).Name );
    for i = 1:length(wvbeg)
      fprintf( 1, ' %d-%d', wvbeg(i), wvend(i) );
    end
    fprintf( 1, '\n' );
  end
end
