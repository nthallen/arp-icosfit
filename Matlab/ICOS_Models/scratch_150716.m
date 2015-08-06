%%
for i=1:length(res)
  check_params(i,res(i));
end
%%
for i=1:length(res2)
  render_model(res2(i));
  title(sprintf('Result %d/%d', i, length(res2)));
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
