%%
% The basic idea here is:
%   Given R1, R2, RR1 and two other parameters, we should be
%   able to solve for the entire RIM/ICOS configuration.
%   Since R1, R2 and RR1 come in discrete values, it makes sense
%   to fix those. For the two additional values to fix, I think
%   it makes sense to fix r1 and Rw1. Generally configurations will
%   scale with r1, so this just establishes the scale of the solution.
%   Rw1 defines the minimum beam spacing to make the Herriott mirror
%   work. Since RL = Rw1/Rs2, it's clear that for a given ICOS
%   configuration, a larger Rw1 will give us a longer RIM, and we'd
%   like to keep it at a minimum.
%% Radius of curvature for plano-convex optics
% These are some of the possible values for R2
IM2 = ispcatalog;
IR2 = unique([IM2.R_cm]);

HM1 = ed_rim_catalog;
HR1 = unique([HM1.R_cm]);
nIR2 = length(IR2);
nHR1 = length(HR1);
ntrials = nIR2*nHR1;

clear P res
P.R1 = 75;
P.r1 = 3;
P.Rw1 = 0.25;
res(ntrials) = exparam(P); % should be empty, just to initialize
trialn = 1;
for IR2i = 1:nIR2
  P.R2 = IR2(IR2i);
  for HR1i = 1:nHR1
    P.RR1 = HR1(HR1i);
    res(trialn) = exparam(P);
    trialn = trialn+1;
  end
end
%%
res = split_results(res);
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
  RD1 = min(RD1(Rr1+0.3 < RD1*2.54/2));
  res(i).RD1 = RD1;
  if isempty(RD1)
    fprintf(1,'Discarding %d: Rr1 exceeds RD1/2\n', i);
    issane(i) = 0;
  end
  
  
  r2 = res(i).r2;
  R2 = res(i).R2;
  D2 = M2Din(M2Rcm == R2);
  D2 = min(D2(r2+0.3 < D2*2.54/2));
  res(i).D2 = D2;
  if isempty(D2)
    fprintf(1,'Discarding %d: r2 exceeds D2/2\n', i);
    issane(i) = 0;
  end
  
  % Radius of curvature cannot exceed mirror radius
  L = res(i).L;
  if L < 5
    fprintf(1,'Discarding %d: cell too short: %f\n', i, L);
    issane(i) = 0;
  end
end
res = res(issane > 0);
%%
%for i = 1:length(res)
% Need to propagate CT1 and CT2 through the catalog
i = 0;
%%
while i < length(res)
  %%
  i = i+1;
  P = ICOS_Model6.props;
  P.R1 = res(i).R1;
  P.R2 = res(i).R2;
  P.r1 = 3*2.54/2;
  P.r2 = res(i).D2*2.54/2;
  P.mirror_spacing = res(i).L;
  P.y0 = res(i).Rr2;
  P.dy = -res(i).Rd2;
  P.dz = res(i).Rs2;
  % P.CT2 = res(i).CT2;
  P.HRC = res(i).RR1;
  P.Hr = res(i).RD1*2.54/2;
  P.HCT = 0.4;
  P.herriott_spacing = res(i).RL;

  P.stop_ICOS = 0;
  P.visible = 1;
  P.visibility = [];
  P.focus = 0;
  P.ICOS_passes_per_injection = 100;
  P.max_rays = 3000;
  P.injection_scale = 1;
  % PM = ICOS_Model6(P);
%% fine tune dy/dz
  P.stop_ICOS = 0;
  P.visible = 0;
  P.visibility = 0;
  P.focus = 0;
  P.HR = 0;
  P.ICOS_passes_per_injection = 100;
  P.max_rays = 3000;
  P.injection_scale = 1;
  P.plot_endpoints = 0;
  P.evaluate_endpoints = 3;
  delta = 2.5e-3;
  best_eccentricity = 2;
  eccentricity = 1;
  iteration = 1;
  f1 = [];
  f2 = [];
  f3 = [];
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
    PM = ICOS_Model6(P,'dy',P.dy+linspace(-delta,delta,11),'dz',P.dz+linspace(-delta,delta,11));
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
        delta = delta/10;
      end
    end
    iteration = iteration + 1;
  end
  %%
  %PM.plot_results('total_power');
  %title(sprintf('Result %d', i));
  %% Now let's see if we can't get the RIM working
  P.stop_ICOS = 1;
  P.HR = 1;
  P.visible = 0;
  P.plot_endpoints = 0;
  P.evaluate_endpoints = 1;
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
    iteration = iteration + 1;
    PM = ICOS_Model6(P,'herriott_spacing', ...
      P.herriott_spacing + linspace(-delta,delta,21));
    %%
    PM.plot_results('total_power');
    title(sprintf('%d/%d: RIM Power: delta = %g iteration %d', ...
      i, length(res), delta, iteration));
    drawnow; shg;
    %%
    [P.herriott_spacing,~,max_power] = PM.identify_minimum('-total_power');
    max_power = abs(max_power);
    delta = delta/10;
  end
  %%
  P.stop_ICOS = 0;
  P.ICOS_passes_per_injection = 20;
  P.max_rays = ceil(max_power+5)*5*P.ICOS_passes_per_injection;
  %%
  P.plot_endpoints = 0;
  P.evaluate_endpoints = 3;
  delta = 1;
  PM = ICOS_Model6(P,'herriott_spacing', ...
    P.herriott_spacing + linspace(-delta,delta,11));
  PM.plot_results('total_power');
  title(sprintf('%d/%d: RIM+ICOS Total Power', i, length(res)));
  drawnow; shg;
  %%
  [P.herriott_spacing,~,~] = PM.identify_minimum('-total_power');
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
