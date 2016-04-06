%%
% HCL_search_B
% This is very much like HCl_search, but without most of the exploration
% code. It just runs straight through.
%
% The current adaptation is to look at whether there are shorter configs
% that we were excluding arbitrarily.
% mnc starts with 'HCl_B4_rd...'
rd = 0.08; % Define the target radius at the detector
th = 13; % Define the target angle from normal at the detector
L = [10 100]; % Cavity length
R2_lim = [-1000 0];

R1 = [25 50 75 100 200 300 400 500 750 1000];
IM2 = ispcatalog;
IR2 = unique([IM2.R_cm]);
if ~isempty(R2_lim)
  IR2 = IR2(IR2 > R2_lim(1) & IR2 < R2_lim(2));
end
R2 = [IR2 R1];
HM1 = ed_rim_catalog;
RR1 = unique([HM1.R_cm]);
%
B = 0.2;
ICOS_margin = B;
mnc = sprintf('HCl_TB%d_rd%02d_th%.1f', round(B*10), floor(rd*100), th);
SR = ICOS_sr_search('mnc',mnc,'L',L,'B',B,'Rw1',0.2,'R1',R1,'R2',R2,...
  'RR1',RR1,'rd',rd,'th',th, 'ICOS_margin', ICOS_margin, 'RD1_margin', 2.0);
SR.enumerate;
%%
SR.design_tolerance;
SR.build_tolerance;
%%
ddBL = diff([SR.Summary.dBL])';
r_max = [SR.Summary.r_max]';
mirror_resolution = 2.54/2;
mirror_margin = SR.SRopt.ICOS_margin;
mirror_lip = 2.54/8; % 1/8"
r_mirror = ceil((r_max+mirror_margin+mirror_lip)/mirror_resolution) * ...
  mirror_resolution - mirror_lip;
% r_mirror = r_max+mirror_margin+mirror_lip; % This is theoretical min
L = [SR.Summary.L]';
Volume = pi*L.*(r_mirror-mirror_lip).^2*1e-3; % liters
%%
% Now pick configurations based on ddBL, V
ddBLok = ddBL > 0.2;
N_configs = 30;
[~,VI] = sort(Volume);
Vlow = sum(cumsum(ddBLok(VI))<=N_configs);
SR.select(VI(1:Vlow));
SR.deselect(find(~ddBLok));
% Skip ahead to savefile

%%
SR.savefile;
%%
SR.focus('focus_visible',0,'max_lens_radius',10);
%%
pattern = ['IS_' SR.SRopt.mnc '*.mat'];
results = dir(pattern);
focuses = zeros(length(results),1);
configs = zeros(length(results),1);
for i=1:length(results)
  load(results(i).name);
  configs(i) = length(IS.res1);
  focuses(i) = length(IS.res2);
end
fprintf(1,'%d possible configurations\n', length(results));
fprintf(1,'%d passed search_ICOS_RIM\n', sum(configs));
fprintf(1,'%d can be focused with 2 or fewer off-the-shelf lenses\n', ...
  sum(focuses>0));
%%
% Discard the IS files with no configurations
v = find(configs==0);
if any(v)
  for vi = v'
    fprintf(1,'Deleting %s\n', results(vi).name);
    delete(results(vi).name);
  end
  results = dir(pattern);
end
%%
% Automatically choose focuses
results = dir(pattern);
Max_Foci = 2;
for i=1:length(results)
  load(results(i).name);
  Ltot = [IS.res2.Ltot];
  Lok = find(Ltot <= min(Ltot)+2.0);
  rpos = zeros(length(Ltot),1);
  for j=1:length(rpos)
    lp = IS.res2(j).Lens_Space;
    rpos(j) = lp(end)/(lp(end)+IS.res2(j).detector_spacing);
  end
  rdist = abs(rpos-0.5);
  [~,RI] = sort(rdist(Lok));
  N = min(length(RI),Max_Foci);
  for j=1:N
    IS.res2(Lok(RI(j))).sel = 1;
  end
  % IS.explore_focus;
  IS.savefile;
end
%%
results = dir(pattern);
for i=1:length(results)
  load(results(i).name);
  IS.analyze('RD1_margin',2);
  close all;
end
%%
res = collect_results('files',pattern,'Rr1',[], ...
  'r2',[],'D1',[],'Lens_Space',[],'detector_spacing',[]);
%%
for i=1:length(res)
  mnc = res(i).mnc;
  % SR#, R1, R2, RR1, m, k
  A = sscanf(mnc,[SR.SRopt.mnc '.%d_%d_%d_%d_%d.%d']);
  res(i).m = A(5);
  res(i).k = A(6);
  r1 = [res.r1];
  L = [res.L];
  res(i).Volume = pi*L(i).*(res(i).D1*2.54/2).^2*1e-3; % liters
  res(i).L_tot = res(i).RL + res(i).L + sum(res(i).Lens_Space) + ...
    res(i).detector_spacing;
  % res(i).Rth = asind(res(i).Rw1/res(i).Rr1);
end
