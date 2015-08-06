%%
% The basic idea here is:
%   Given R1, R2, RR1 and two other parameters, we should be
%   able to solve for the entire RIM/ICOS configuration.
%   Since R1, R2 and RR1 come in discrete values, it makes sense
%   to fix those. For the two additional values to fix, I thought
%   it would makes sense to fix r1 and Rw1, but that didn't really
%   work out so well because it produced cell lengths L that differed
%   significantly from what we have available. The plan here is to
%   fix R1, R1, RR1, L and Rw1.
%
%   After calculating theoretical parameters via exparam, This function
%   tests those parameters via ray tracing and optimizes them to correct
%   for some of the approximations used in the simplified theoretical
%   model.
%% Radius of curvature for plano-convex optics
% These are some of the possible values for R2
IM2 = ispcatalog;
IR2 = unique([IM2.R_cm]);
IR2 = 75;

HM1 = ed_rim_catalog;
HR1 = unique([HM1.R_cm]);
nIR2 = length(IR2);
nHR1 = length(HR1);
ntrials = nIR2*nHR1;

clear P res
P.R1 = 75;
P.L = 50; %3;
P.Rw1 = 0.25;
Res = cell(ntrials,0);
% res(ntrials) = exparam1(P); % should be empty, just to initialize
trialn = 1;
for IR2i = 1:nIR2
  P.R2 = IR2(IR2i);
  for HR1i = 1:nHR1
    P.RR1 = HR1(HR1i);
    Res{trialn} = exparam(P)';
    trialn = trialn+1;
  end
end
%%
res = [Res{:}]';
%%
nres = length(res);
issane = ones(nres,1);
HM1Rcm = [HM1.R_cm];
HM1Din = [HM1.dia_in];
M2Rcm = [IM2.R_cm];
M2Din = [IM2.dia_in];
for i = 1:nres
  % Now let's check for some sanity:
  % Is Rr1+0.3 < RD1? (spot radius less than mirror radius)
  % is RR1 < RD1? (radius of curvature less than mirror radius)
  Rr1 = res(i).Rr1;
  RR1 = res(i).RR1;
  RD1 = HM1Din(HM1Rcm == RR1);
  RD1a = min(RD1(Rr1+0.3 < RD1*2.54/2));
  res(i).RD1 = RD1a;
  if isempty(RD1a)
    fprintf(1,'Discarding %d: Rr1 (%.1f) exceeds RD1/2 (%.1f)\n', i, Rr1, max(RD1)*2.54/2);
    issane(i) = 0;
  end
  
  r2 = res(i).r2;
  R2 = res(i).R2;
  D2 = M2Din(M2Rcm == R2);
  D2a = min(D2(r2+0.3 < D2*2.54/2));
  res(i).D2 = D2a;
  if isempty(D2)
    fprintf(1,'Discarding %d: r2 (%.1f) exceeds D2/2 (%.1f)\n', i, r2, max(D2)*2.54/2);
    issane(i) = 0;
  end
  
  r1 = res(i).r1;
  D1 = 3; % fixed for now
  if r1 > D1*2.54/2
    fprintf(1,'Discarding %d: r1 (%.1f) exceeds D1/2 (%.1f)\n', i, r1, D1*2.54/2);
    issane(i) = 0;
  end
  
  L = res(i).L;
  if L < 5
    fprintf(1,'Discarding %d: cell too short: %f\n', i, L);
    issane(i) = 0;
  end
end
%%
res = res(issane > 0);
%%
%for i = 1:length(res)
% Need to propagate CT1 and CT2 through the catalog
f1 = [];
f2 = [];
f3 = [];
i = 0;
%%
while i < length(res)
  %%
  i = i+1;
  P = render_model(res(i));
%% fine tune dy/dz
  P.stop_ICOS = 0;
  P.visible = 0;
  P.visibility = 0;
  P.focus = 0;
  P.HR = 0;
  P.ICOS_passes_per_injection = 100;
  P.max_rays = 3000;
  P.plot_endpoints = 0;
  P.evaluate_endpoints = 3;
  P.skip.overlap = 1;
  P.skip.total_power = 1;
  P.skip.mean_angle = 1;
  P.skip.RIM_passes = 1;
  delta = 2.5e-3;
  best_eccentricity = 2;
  eccentricity = 1;
  iteration = 1;
  %%
  if isempty(f1)
    f1 = figure;
    pos = get(f1,'position');
    pos(1) = 0;
    set(f1,'position',pos);
  else
    figure(f1);
  end
  while delta > 1e-4 && eccentricity < best_eccentricity && iteration < 10
    %%
    best_eccentricity = eccentricity;
    PM = ICOS_Model6(P,'dy',P.dy+linspace(-delta,delta,5),'dz',P.dz+linspace(-delta,delta,5));
    PM = PM.clean_results;
    PM.plot_results('eccentricity');
    title(sprintf('%d/%d: ICOS optimization delta = %g, iteration %d', ...
      i, length(res), delta, iteration));
    drawnow; shg;
    [new_dy,new_dz,eccentricity,ri,ci] = PM.identify_minimum('eccentricity');
    if isempty(new_dy) || isempty(new_dz)
      delta = delta*10;
    elseif eccentricity < best_eccentricity
      P.dy = new_dy;
      P.dz = new_dz;
      if ri ~= 1 && ri ~= size(PM.Results.eccentricity,1) && ...
        ci ~= 1 && ci ~= size(PM.Results.eccentricity,2)
        delta = delta/2;
      end
    end
    iteration = iteration + 1;
  end
  %% Now let's see if we can't get the RIM working
  P.stop_ICOS = 1;
  P.HR = 1;
  P.visible = 0;
  P.plot_endpoints = 0;
  P.evaluate_endpoints = 1;
  P.skip.RIM_passes = 0;
  delta = 2; % cm +/-
  if isempty(f2)
    f2 = figure;
    pos(1) = pos(1) + pos(3);
    set(f2,'position',pos);
  else
    figure(f2);
  end
  iteration = 0;
  %%
  while delta > 0.1
    %%
    iteration = iteration + 1;
    PM = ICOS_Model6(P,'herriott_spacing', ...
      P.herriott_spacing + linspace(-delta,delta,21));
    %%
    criteria = 'eccentricity';
    PM = PM.clean_results;
    PM.Results.eccentricity(PM.Results.RIM_passes <= 1) = NaN;
    PM.plot_results(criteria);
    title(sprintf('%d/%d: RIM eccentricity: delta = %g iteration %d', ...
      i, length(res), delta, iteration));
    drawnow; shg;
    %%
    [P.herriott_spacing,~,RIM_passes] = PM.identify_minimum('-RIM_passes');
    RIM_passes = abs(RIM_passes)+1;
    [P.herriott_spacing,~,~] = PM.identify_minimum(criteria);
    delta = delta/10;
  end
  %%
  P.stop_ICOS = 0;
  P.ICOS_passes_per_injection = 20;
  P.max_rays = ceil(RIM_passes+5)*5*P.ICOS_passes_per_injection;
  %%
  P.plot_endpoints = 0;
  P.evaluate_endpoints = 3;
  delta = .1;
  PM = ICOS_Model6(P,'herriott_spacing', ...
    P.herriott_spacing + linspace(-delta,delta,11));
  PM = PM.clean_results;
  PM.Results.eccentricity(PM.Results.RIM_passes <= 1) = NaN;
  PM.plot_results(criteria);
  title(sprintf('%d/%d: RIM+ICOS %s', i, length(res), criteria));
  drawnow; shg;
  [P.herriott_spacing,~,~] = PM.identify_minimum(criteria);
  %%
  P.visible = 1;
  P.visibility = [];
  P.focus = 1;
  P.evaluate_endpoints = -1;
  if isempty(f3)
    f3 = figure;
    pos(1) = pos(1) + pos(3);
    set(f3,'position',pos);
  else
    figure(f3);
  end
  PM = ICOS_Model6(P);
  title(sprintf('Result %d', i));
  drawnow; shg;
  %%
  res(i).ORd2 = -P.dy;
  res(i).ORs2 = P.dz;
  res(i).ORL = P.herriott_spacing;
end
%%
save exexparam3a_75.mat res