function varargout = get_waveform_params( waveform, varargin )
% [ val1, val2 ] = get_waveform_params( waveform, var1, var1_default, var2, ...
%                                var2_default, ... )
% vals = get_waveform_params( waveform, var1, var1_def, ...
%                              var2, var2_def, ... )
% vals = get_waveform_params( waveform, 'struct', var1, var1_def )
%
% Defaults for at least one parameter must be specified.

% Parameters currently in use include:
% matchline6a:
%   'hlf', 24
%   'mfhw', 1
%   'cv_hlf', 10
%   'cv_mfhw', 1
%   'holdoff', 4e-4
%   'SignalRegion', [round(wv.TzSamples + holdoff*SampleRate): ...
%                          wv.NetSamples - wv.TzSamples - 1]
% etln_fit7 & etln_fit8:
%   'prefilterwidth', 5
%   'X', [10.2817    48.3    0   -2.6921  .1645924   -3.7796   .0689779 ]
%   'SignalRegion', [wv.TzSamples+100:wv.NetSamples-wv.TzSamples-20]
%   'threshold', .07
%  (etln_fit7 used to use 'range' instead of 'SignalRegion')
%
% sensitivity:
%   'holdoff', 4e-4
%   'SignalRegion', []
%
% waves_editor:
%   'RampRegions'
%   'CurRamp'

argi = 1;
if nargout == 1 && length(varargin) >= 3
  structout = 1;
  if strcmp(varargin{1},'struct')
    argi = 2;
  end
else
  if length(varargin) ~= nargout*2
    error('Number of output arguments must match inputs plus defaults');
  end
  structout = 0;
end
fname = findinpath( [ waveform '_etln.mat' ], { '.', '..', '../..' } );
if length(fname)
  fprintf(1, 'Reading waveform configuration from %s\n', fname );
  vals = load(fname);
else
  vals = struct( varargin{argi}, varargin{argi+1} );
end
if structout
  for i=argi+1:2:length(varargin)
    fld = varargin{i-1};
    if ~isfield(vals,fld)
      vals.(fld) = varargin{i};
    end
  end
  varargout{1} = vals;
else
  for i = 1:nargout
    fld = varargin{2*i-1};
    if isfield(vals,fld)
      varargout{i} = vals.(fld);
    else
      varargout{i} = varargin{2*i};
    end
  end
end