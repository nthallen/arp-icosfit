function Analyze_IM6b(res, prefix, varargin)
% Analyze_IM6b(res, prefix [, options] );
%   res: output from exexparam4
%   prefix: string to prefix output files
%   select: optional vector of indices into res
nres = length(res);
Opt.select = 1:nres;
Opt.Npasses = 100;
Opt.Nsamples = 100;
for i=1:2:length(varargin)-1
  if isfield(Opt, varargin{i})
    Opt.(varargin{i}) = varargin{i+1};
  else
    error('MATLAB:HUARP:InvalidOption', ...
      'Invalid option: "%s"', varargin{i});
  end
end
if iscolumn(Opt.select)
  Opt.select = Opt.select';
end
ff = [];
for i=Opt.select
  ofile = sprintf('%s.%d_%dx%d.mat', prefix, i,Opt.Npasses,Opt.Nsamples);
  if ~exist(ofile, 'file')
    P = render_model(res(i));
    P.visible = 0;
    % P.HR = 0;
    P.focus = 1;
    P.evaluate_endpoints = -1;
    opt_n = 4 + length(P.Lenses);
    IB = ICOS_beam(@ICOS_Model6, P);
    IB.Sample('beam_samples', Opt.Nsamples, ...
      'ICOS_passes', Opt.Npasses, 'opt_n', opt_n, ...
      'Track_Power', 1);
    ff2 = IB.Integrate;
    if ~isempty(ff)
      delete(ff);
    end
    ff = ff2;
    save(ofile, 'IB');
    fprintf(1, 'Saved result to %s\n', ofile);
  else
    fprintf(1, 'Skipping result %d\n', i);
  end
end
