%%
% HTW_search.m
% Need to accomodate possible change in n for 3765-3780 cm-1
%   Somewhere in the range of 2.4465 to 2.4467 based on interpolation
rd = 0.08;
th = 13;
B = 0.2;
ICOS_margin = 0.2;
mnc = 'HTW'; % sprintf('HTW_rd%02d_th%.1f', floor(rd*100), th);
SR = ICOS_sr_search('mnc',mnc,'B',B,'Rw1',0.15, 'rd', rd, 'th', th, ...
  'R1',[75 100 150],'R2',[75 100 150],'L',[40 100],'optics_n',2.4466, ...
  'ICOS_margin', ICOS_margin, 'RD1_margin', 2.0);
SR.enumerate;
%%
SR.design_tolerance;
SR.build_tolerance;
%%
ddBL = diff([SR.Summary.dBL])';
r_max = [SR.Summary.r_max]';
mirror_resolution = 2.54/2;
mirror_margin = ICOS_margin;
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
%%
SR.explore_build_tolerance; % BT vs RLmin
SR.explore_build_tolerance(2); % BT vs DT
SR.explore_build_tolerance(3); % BT vs r2/r1
% SR.explore_build_tolerance(4); % R1 vs R2 (uninteresting)
%
SR.explore_build_tolerance(5); % L
%%
SR.savefile;
%%
L = [SR.Summary.L]';
SR.deselect(1:length(L));
SR.select(find(abs(L-56.2849) < .001));
%%
SR.focus('focus_visible',0);
%%
% Here we discard most of what we just did, but reuse the IS.res1 solution:
IS.res2 = [];
% These are the lenses we have for HTW:
fix = {'ZC_PM_50_50','ZC_PM_25_50' };
th = 13; % 13 degrees
IS.search_focus2('max_lenses',2,'fix_lenses',fix,'det_acc_limit',th,...
  'allow_nondecreasing_focus',1);
%%
IS.explore_focus;
IS.savefile;
%%
results = dir('IS_HTW*.mat');
for i=1:length(results)
  load(results(i).name);
  fprintf(1,'%d: %d focus solutions: %s \n', i, length(IS.res2), results(i).name);
  fprintf(1,'%d: d2 = %.4f\n', i, IS.res1(1).d2);
  check_params(i, IS.res1(1));
  for j=1:length(IS.res2)
    fprintf(1, '  %d: ', j);
    Lenses = IS.res2(j).Lenses;
    for k=1:length(Lenses)
      fprintf(1, ' %s', Lenses{k});
    end
    fprintf(1,'\n');
  end
end
%%
results = dir('IS_HTW*.mat');
for i=1:length(results)
  load(results(i).name);
  fprintf(1,'%d: RL: %.1f R1: %.0f R2: %.0f L: %.1f r2/r1: %.1f\n', i, ...
    IS.res1.RL, IS.res1.R1, IS.res1.R2, IS.res1.L, IS.res1.r2/IS.res1.r1);
end
%%
% Automatically choose focuses
pattern = ['IS_' SR.SRopt.mnc '*.mat'];
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
pattern = ['IS_' SR.SRopt.mnc '*.mat'];
results = dir(pattern);
for i=1:length(results)
  load(results(i).name);
  IS.analyze('RD1_margin',2);
  close all;
end
%%
% This assumes IB has been loaded
IB.draw('r1',2.54,'r2',2.54,'HR',0,'CT1',0.5,'CT2',0.5,'visibility', 0);
title('HTW ICOS Configuration: no RIM');
view(0,90);
%%
xlabel('cm');
ylabel('cm');
%%
fname = [IB.IBP.mnc '.png'];
print(gcf, '-dpng', fname);
%%
PM = IB.draw('r1',2.54,'r2',2.54,'HR',0,'CT1',0.5,'CT2',0.5,'visibility', [0 0 0]);
title('HTW ICOS Focus: no RIM');
view(0,90);
%%
xlabel('cm');
ylabel('cm');
%%
fname = [IB.IBP.mnc '_focus.png'];
print(gcf, '-dpng', fname);
%%
fprintf(1,'\nOptical Configuration Summary for %s\n\n', IB.IBP.mnc);
X = IB.P.mirror_spacing;
fprintf(1,'Space between M1 and M2 (face to face): %.2f cm\n',X);
X = X + IB.P.CT2+IB.P.Lens_Space(1);
fprintf(1,'L1 leading position: %.2f cm\n', X);
fprintf(1,'Space between M2 and L1: %.2f cm\n',IB.P.Lens_Space(1));
if length(IB.P.Lens_Space) == 2
  X = X + PM.M.Optic{4}.CT+IB.P.Lens_Space(2);
  fprintf(1,'L2 leading position: %.2f cm\n', X);
  fprintf(1,'Space between L1 and L2: %.2f cm\n',IB.P.Lens_Space(2));
  X = X + PM.M.Optic{5}.CT+IB.P.detector_spacing;
else
  X = X + PM.M.Optic{4}.CT+IB.P.detector_spacing;
end
fprintf(1,'Detector position: %.2f cm\n', X);
fprintf(1,'Space between L%d and detector: %.2f cm\n', ...
  length(IB.P.Lens_Space), IB.P.detector_spacing);
%%
P = IB.P;
P.r1 = 2.54; % 2" diamter
P.r2 = 2.54; % 2" diameter
P.Hr = 2.54; % 2" diameter: what the heck
% This is target for mirror in free mount
%draw_targets(P,[4*2.54, 0.345]);
% This is target for RIM + Carolina's cavity end cap
draw_targets(P,[IB.P.herriott_spacing, 2.78]);
%%
fname = [IB.IBP.mnc '_target.png'];
print(gcf,'-dpng',fname,'-r600');
