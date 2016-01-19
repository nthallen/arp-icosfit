%% scratch_160113: Model of benchtop HCl/ME system with crappy mirrors
R = 75;
L = 50;
w_r = sqrt(L*(2*R-L)/R^2);
phi = asind(w_r);
r1 = 1*2.54;
w1 = r1 * w_r;
B = 0.4;
Rw1 = 1;
s1 = w1/L;
n = 2.4361;
Rs2 = s1;
RL = Rw1/Rs2;
SP = struct('R1',R,'R2',R,'L',L,'RL',RL,'Rw1',Rw1);
Res = exparam(SP);
% IS = ICOS_search('R1',R,'R2',R,'L',L,
IS = ICOS_search('mnc', 'JW','R1',SP.R1,'R2',SP.R2,'L',SP.L, ...
  'RR1',Res.RR1,'Rw1',SP.Rw1, 'RL_lim', [0.90,1.1]*SP.RL);
%%
IS.search_ICOS_RIM;
%%
IS.search_focus2('max_lenses',2);
%%
V = find([IS.res2.sel]);
ff = [];
for ni = 1:length(V)
  i = V(ni);
  P = render_model(IS.res2(i));
  P.HR = 0; % No RIM for now
  P.T = .005; % 5000 ppm
  Opt.ICOS_passes = 50;
  Opt.Nsamples = 100;
  IBopt = {};
  
  n_optics = 4 + length(P.Lenses);
  opt_n = n_optics;
  Track_Power = 1;
  IBmnc = sprintf('%s.%d_%dx%d', IS.ISopt.mnc, ...
    IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
  ofile = sprintf('IB_%s.mat', IBmnc);
  if ~exist(ofile, 'file')
    P.visible = 0;
    % P.HR = 0;
    P.focus = 1;
    P.evaluate_endpoints = -1;
    IB = ICOS_beam(@ICOS_Model6, P);
    IB.Sample('beam_samples', Opt.Nsamples, ...
      'ICOS_passes', Opt.ICOS_passes, 'opt_n', opt_n, ...
      'n_optics', n_optics, 'Track_Power', Track_Power, ...
      'mnc', IBmnc, IBopt{:});
    if opt_n == n_optics
      ff2 = IB.Integrate;
      if ~isempty(ff)
        delete(ff);
        drawnow;
      end
      ff = ff2;
    end
    IB.savefile;
    % save(ofile, 'IB');
    % fprintf(1, 'ICOS_beam %d Saved result to %s\n', i, ofile);
    Pwr = IB.PowerSummary;
    IS.res2(i).NH = Pwr.NH;
    if isfield(Pwr,'max_pwr')
      IS.res2(i).max_pwr = Pwr.max_pwr;
    elseif ~isfield(IS.res2,'max_pwr')
      IS.res2(i).max_pwr = [];
    end
    IS.savefile;
  else
    fprintf(1, 'Skipping result %d\n', i);
  end
end
