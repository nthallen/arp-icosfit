function [line_pos, Sout, wvout] = get_line_position(region, suffix, PTEfile)
% [line_pos, S, wv] = get_line_position(region, suffix, PTEfile)
% Inputs:
%  region, suffix: For specifying base ICOSout.region.suffix
%  PTEfile: defaults to S.PTEfile
% Outputs:
%  linepos: sample number of line centers
%  S: Output structure of ICOS_setup
%  wv: Waveform structure
S = ICOS_setup(['ICOSout.' region '.' suffix ]);
if nargin < 3 || isempty(PTEfile)
    PTEfile = S.PTEfile;
end
scans = S.scannum;
wv = waves_used(scans);
if length(wv) > 1
    error('More than one waveform in scan region');
end
SignalRegion = get_waveform_params( wv.Name, ...
    'SignalRegion', wv.TzSamples+100:wv.NetSamples-wv.TzSamples-20 );
nu = get_nu(region, suffix, SignalRegion, PTEfile);
line_pos = zeros(size(S.nu_P));
for i=1:length(scans)
    line_pos(i,:) = interp1(nu(:,i),SignalRegion,S.nu_P(i,:));
end
if nargout > 1
    Sout = S;
    if nargout > 2
        wvout = wv;
    end
end
