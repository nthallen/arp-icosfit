%%
for i=1:length(res)
  check_params(i,res(i));
end
%%
for i=1:length(res)
  render_model(res(i));
  title(sprintf('Result %d/%d', i, length(res)));
  drawnow; shg;
end
%%
res = exexparam4('exexparam3_2save.mat');
%%
Analyze_IM6b(res,'AIM6btmp.2', 'select', 2);
%%
Analyze_IM6b(res,'AIM6btmp.2', 'select', 2, 'Nsamples', 500, 'ICOS_passes', 10);
%%
Analyze_IM6b(res,'AIM6btmp2.2', 'select', 2, 'Nsamples', 100, ...
  'ICOS_passes', 1, 'rng_state', IB.IBP.rng_state, 'opt_n', 1, 'Npass_dim', 1);
%%
Analyze_IM6b(res,'AIM6btmp3.2', 'select', 2, 'Nsamples', 100, ...
  'ICOS_passes', 1, 'rng_state', IB.IBP.rng_state, 'opt_n', 1, 'Npass_dim', 1, 'Track_Power', 1);
%%
% AIM6b.2 is based on exexparam3/4 with r1 = 2 cm.
Analyze_IM6b(res, 'AIM6b.2', 'Npasses', 10000,'select',1:50);
%%
res2 = exexparam4('exexparam3a_75.mat');
%%
% AIM6b.3a.75 is from exexparam3a/4 with L = 50cm Rw1 = 0.25 and R2=75cm
% (The selection here is limited to Ltot < 1m)
Analyze_IM6b(res2, 'AIM6b.3a.75', 'select', 1:18, 'Track_Power', 1);
%%
% This is from exexparam3a_75_5 => exexparam4_75_5a
Analyze_IM6b(res, 'AIM6b.3a.75_5a', 'select', 8, 'Track_Power', 1);
%%
% This is from exexparam3a_w5_L50 => exexparam4_w5_L50
Analyze_IM6b(res, 'AIM6b_w75_L50', 'select', 9, 'Track_Power', 1);
%%
% This is from exexparam3a_w25_L50 => exexparam4_w25_L50
Analyze_IM6b(res, 'AIM6b_w25_L50', 'select', 32, 'Track_Power', 1);
%%
% Start testing of ICOS_search
IS = ICOS_search('R1', 75, 'Rw1', 0.4, 'mnc', 'w4_L50', 'L', 50, 'RL_lim', [-inf 50]);
IS.search_ICOS_RIM;
IS.search_focus('select1', 3);
IS.analyze('select2',3);
%%
IS = ICOS_search('R1', 75, 'Rw1', 0.25, 'mnc', 'w25_L50', 'L', 50, 'RL_lim', [-inf 50]);
IS.search_ICOS_RIM;
%%
IS.search_focus('select', 5);
%%
IS.analyze('select',[2:5]);
%%
% This is specifically looking at M2=ZC-PX-38-200, R2=-28.12
IS = ICOS_search('R1', 75, 'Rw1', 0.25, 'mnc', 'w25_L50_R28', 'L', 50, ...
  'R2_lim', [-29 -28], 'RL_lim', [2 50]);
IS.search_ICOS_RIM;
%%
IS.search_focus('select', 1);
%%
IS.analyze('select',[1:2]);
%%
% 75/75 w30
IS = ICOS_search('R1', 75, 'Rw1', 0.3, 'mnc', 'w30_L50', 'L', 50, ...
  'RL_lim', [2 50]);
IS.search_ICOS_RIM;
%%
IS.search_focus('select', [6 8]);
%%
IS.analyze('select',[3 8]);

%% Play with exparam to see how w1 affects L, RL
nonempty = {};
for Rw1 = 26:29
  mnc = sprintf('ex_R28_w%d_r24', Rw1); % 'ex_' prefix to indicate experimental
  IS = ICOS_search('R1', 75, 'Rw1', Rw1/100, 'mnc', mnc, ...
    'R2_lim', [-29 -28], 'r1', 2.4, ...
    'RL_lim', [2 100]);
  IS.search_ICOS_RIM;
  if ~isempty(IS.res1)
    nonempty{end+1,1} = mnc;
  end
end
nonempty
%%
nonempty = {};
for Rw1 = 26:29
  ifile = sprintf('IS_R28_w%d_r24.mat', Rw1); % 'ex_' prefix to indicate experimental
  load(ifile);
  IS.search_focus;
  if ~isempty(IS.res2)
    nonempty{end+1,1} = mnc;
  end
end
nonempty
%%
for i=1:190
  ofile = sprintf('AIM6b.2.%d_10000x100.mat', i);
  if exist(ofile, 'file')
    load(ofile);
    ff = IB.Integrate;
    pause;
    delete(ff);
  end
end
%%
Pwr = []; % This needs to be done once before any processing
% The following then needs to happen after each beam is processed
Ri = 2:PM.M.n_rays;
pre_opt = [PM.M.Rays([PM.M.Rays(Ri).n_inc]).n_opt];
cur_opt = [PM.M.Rays(Ri).n_opt];
A = [pre_opt' cur_opt'];
[B,ia,ib] = unique(A,'rows');
%%
RP = zeros(size(Ri));
inside = zeros(size(Ri'));
for i=1:length(Ri)
  ray = PM.M.Rays(Ri(i)).ray;
  RP(i) = ray.P;
  inside(i) = ray.Inside;
end
outside = ~inside;
%%
for i=1:length(B)
  fld = sprintf('R%d_%d', B(i,1), B(i,2));
  if ~isfield(Pwr, fld)
    Pwr.(fld) = struct('I',0,'O',0);
  end
  for pre = 'IO'
    if pre == 'I'
      icond = inside;
    else
      icond = outside;
    end
    Pwr.(fld).(pre) = Pwr.(fld).(pre) + sum(RP((ib == i) & icond));
  end
end

%%
% Compare results for exexparam3 and exexparam3a
load('exexparam3_2save.mat');
clear P
i = 0;
%%
while i < length(res)
  %%
  i = i+1;
  %%
  P.R1 = res(i).R1;
  P.L = res(i).L;
  P.Rw1 = res(i).Rw1;
  P.R2 = res(i).R2;
  P.RR1 = res(i).RR1;
  res2 = exparam(P);
  %%
  check_solution(res(i));
  %%
  check_params(1, res(i));
  %%
  check_params(2, res2);
end

