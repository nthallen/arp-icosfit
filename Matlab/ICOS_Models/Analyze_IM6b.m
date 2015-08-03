function Analyze_IM6b(res, prefix, varargin)
% Analyze_IM6b(res, prefix [, options] );
%   res: output from exexparam4
%   prefix: string to prefix output files
%   select: optional vector of indices into res
nres = length(res);
Opt.select = 1:nres;
Opt.ICOS_passes = 50;
Opt.Nsamples = 100;
% Opt.Herriott_passes = 100;
% Opt.n_optics = 5;
% Opt.rng_state = [];
Opt.opt_n = [];
IBopt = {};
for i=1:2:length(varargin)-1
  if isfield(Opt, varargin{i})
    Opt.(varargin{i}) = varargin{i+1};
  else
    IBopt{end+1} = varargin{i};
    IBopt{end+1} = varargin{i+1};
  end
end
if iscolumn(Opt.select)
  Opt.select = Opt.select';
end
ff = [];
for i=Opt.select
  ofile = sprintf('%s.%d_%dx%d.mat', prefix, i,Opt.ICOS_passes,Opt.Nsamples);
  if ~exist(ofile, 'file')
    P = render_model(res(i));
    P.visible = 0;
    % P.HR = 0;
    P.focus = 1;
    P.evaluate_endpoints = -1;
    n_optics = 4 + length(P.Lenses);
    if ~isempty(Opt.opt_n)
      opt_n = Opt.opt_n;
    else
      opt_n = n_optics;
    end
    IB = ICOS_beam(@ICOS_Model6, P);
    IB.Sample('beam_samples', Opt.Nsamples, ...
      'ICOS_passes', Opt.ICOS_passes, 'opt_n', opt_n, ...
      'n_optics', n_optics, IBopt{:});
    if opt_n == n_optics
      ff2 = IB.Integrate;
      if ~isempty(ff)
        delete(ff);
      end
      ff = ff2;
    end
    save(ofile, 'IB');
    fprintf(1, 'Saved result to %s\n', ofile);
  else
    fprintf(1, 'Skipping result %d\n', i);
  end
end
