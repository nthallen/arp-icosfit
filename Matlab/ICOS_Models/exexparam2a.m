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

% HM1 = ed_rim_catalog;
% HR1 = unique([HM1.R_cm]);
nIR2 = length(IR2);
r2 = [1:.1:3];
nr2 = length(r2);
ntrials = nIR2*nr2;

clear P res
P.R1 = 75;
P.r1 = 3;
P.Rw1 = 0.25;
P.L = 50;
res(ntrials) = exparam(P); % should be empty, just to initialize
trialn = 1;
for IR2i = 1:nIR2
  P.R2 = IR2(IR2i);
  for r2i = 1:nr2
    P.r2 = r2(r2i);
    resl = exparam(P);
    res(trialn) = resl;
    trialn = trialn+1;
  end
end
%%
res = split_results(res);
%%
nres = length(res);
issane = ones(nres,1);
% HM1Rcm = [HM1.R_cm];
% HM1Din = [HM1.dia_in];
M2Rcm = [IM2.R_cm];
M2Din = [IM2.dia_in];
for i = 1:nres
  % Now let's check for some sanity:
  % Is Rr1+0.3 < RD1? (spot radius less than mirror radius)
  % is RR1 < RD1? (radius of curvature less than mirror radius)
  Rr1 = res(i).Rr1;
  RD1 = (Rr1+0.3)*2/2.54;
  res(i).RD1 = RD1;
  if isempty(RD1)
    issane(i) = 0;
  end
  
  
  r2 = res(i).r2;
  R2 = res(i).R2;
  D2 = M2Din(M2Rcm == R2);
  D2 = min(D2(r2+0.3 < D2*2.54/2));
  res(i).D2 = D2;
  if isempty(D2)
    issane(i) = 0;
  end
  
  % Radius of curvature cannot exceed mirror radius
end
%%
res = res(issane > 0);
%%
% for i = 1:length(res)
  i = 1;
%%
if i < length(res)
  P = ICOS_Model6.props;
  P.R1 = res(i).R1;
  P.R2 = res(i).R2;
  P.r1 = 3*2.54/2; % Mirror radius (not beam spot radius)
  P.r2 = res(i).D2*2.54/2;
  P.mirror_spacing = res(i).L;
  P.y0 = res(i).Rr2;
  P.dy = -res(i).Rd2;
  P.dz = res(i).Rs2;
  % P.CT2 = res(i).CT2;
  P.HRC = res(i).RR1;
  P.Hr = res(i).RD1*2.54/2;
  P.HR = 0;
  % P.HCT = 0.4;
  P.herriott_spacing = res(i).RL;
  P.visible = 1;
  P.visibility = [];
  P.focus = 0;
  P.ICOS_passes_per_injection = 100;
  P.max_rays = 3000;
  P.injection_scale = 1;
  PM = ICOS_Model6(P);
  title(sprintf('Result %d', i));
  i = i+1;
end
  %%
%  pause;
%end
