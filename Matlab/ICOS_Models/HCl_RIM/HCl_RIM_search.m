%%
% HCL_RIM_search
% This is very much like HCl_search, but without most of the exploration
% code. It just runs straight through.
%
% The current adaptation is to look at whether there are shorter configs
% that we were excluding arbitrarily.
% mnc starts with 'HCl_B4_rd...'
rd = 0.08; % Define the target radius at the detector
th = 14; % Define the target angle from normal at the detector
% Cavity length is nominally 50 cm, but there is an additional 1.42 cm
% in each end cap to the front edge of the mirror, and each mirror
% has a sag of 0.10 cm (3" mirror with RoC of 75 cm), so the total
% distance is 50 + 2*(1.42 + 0.10) = 50 + 2*1.52 = 50 + 3.04 = 53.04 cm
% distance is 50 + 2*(0.1 + 0.1) = 50.4 cm
L_lim = [30 90]; % Cavity length
% R2_lim = [-1000 0];
R2_lim = [];

R1 = 200;
% R1 = 75;
% IM2 = ispcatalog;
% IR2 = unique([IM2.R_cm]);
% if ~isempty(R2_lim)
%   IR2 = IR2(IR2 > R2_lim(1) & IR2 < R2_lim(2));
% end
% R2 = [IR2 R1];
R2 = 1000;
% HM1 = ed_rim_catalog;
% RR1 = unique([HM1.R_cm]);
%
RR1 = 100; % What we have in stock
B = 0.2;
Rw1 = B*3/2; %look into further
ICOS_margin = 0.2;
mnc = sprintf('HCl_RIM100_B%d_Rw%1d_rd%02d_th%.1f', round(B*10), ...
  round(Rw1*10), floor(rd*100), th);
%%
% This search is highly constrained based on previous exploratory
% searches. At this point, we have selected the RIM, M1 and M2, so the
% search space is quite small. We have also provided locally modified
% versions of isp_catalog.m and custom_lenses.m to limit the focusing
% search to optics we have already purchased.
SR = ICOS_sr_search('mnc',mnc,'L',L_lim,'R1',R1,'R2',R2,'Rw1',Rw1, ...
  'B', B, 'rd', rd, 'th', th, 'RR1', RR1);
SR.enumerate
%%
SR.build_tolerance
SR.design_tolerance
SR.explore_build_tolerance
%%
% By examination, this is the one that match's Nick's summary
SR.select(14);
%%
SR.focus('focus_visible',0,'max_lens_radius',10);
SR.savefile;
%%
% RD1_margin selected to push RIM to 3". r1 & r2 set for 2" mirrors minus a
% 1/8" margin. Mirror transmission is based on measurements from 6/13/17
IS.analyze('HR',.98,'T',250e-6,'RD1_margin',1.8,'r1',7*2.54/8,'r2',7*2.54/8);
%%
IS = ICOS_search('mnc',mnc,'L_lim',L_lim,'R1',R1,'R2',R2,'Rw1',Rw1, ...
  'beam_diameter', B, 'RL_lim',[5, 50]);
IS.search_ICOS_RIM;
%%
IS.res2 = [];
%%
IS.search_focus2('focus_visible', 1, 'max_lenses', 2, ...
  'allow_nondecreasing_focus',1,'max_focus_length',40, 'dx_min', 1.554);
%%
Ltot = [IS.res2.Ltot];
Max_Foci = 2; % The number of focus solutions to consider
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
IS.explore_focus;
IS.savefile;
%%
% Mirror transmission set high based on previous measurments
IS.analyze('HR',.98,'T',.005,'RD1_margin',2,'r1',3*2.54/2,'r2',3*2.54/2);

%%
pattern = 'IB_HCl_L50_B2*.mat';
res = collect_results('files',pattern,'Rr1',[], ...
  'r2',[],'R1',[],'D1',[],'Lens_Space',[],'detector_spacing',[]);
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

%%
PM = IB.draw('HR',0,'visibility', [0]);
title('HTW ICOS Focus: no RIM');
view(0,90);
%%
fname = [IB.IBP.mnc '.png'];
print(gcf, '-dpng', fname);
%%
PM = IB.draw('HR',0,'visibility', [0 0 0]);
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
fprintf(1,'M1: Diameter: %.2f" RoC: %.1f cm\n', IB.P.r1*2/2.54, IB.P.R1);
fprintf(1,'M2: Diameter: %.2f" RoC: %.1f cm\n', IB.P.r2*2/2.54, IB.P.R2);
for i=1:length(IB.P.Lenses)
  fprintf(1,'L%d: %s\n', i, IB.P.Lenses{i});
end
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
P.r1 = 1.5*2.54; % 3" diamter
P.r2 = 1.5*2.54; % 3" diameter
P.Hr = 1.5*2.54; % 2" diameter: what the heck
% This is target for mirror in free mount
%draw_targets(P,[4*2.54, 0.345]);
% This is target for RIM + Carolina's cavity end cap
% draw_targets(P,[6.00, 0.32]);
% This is target for RIM + HTW end cap
draw_targets(P,[6.00, 1.554]);
%%
fname = [IB.IBP.mnc '_target.png'];
print(gcf,'-dpng',fname,'-r600');

%%
% Redo ICOS_beam analysis with beam_diameter and beam_divergence as
% predicted by the telescope model
IS.ISopt.mnc = 'HCl_L50_B16_Rw3_rd08_th14.0_D12';
IS.ISopt.beam_diameter = 1.6;
IS.analyze('HR',0,'T',.005,'RD1_margin',2,'r1',3*2.54/2,'r2',3*2.54/2, ...
  'beam_diameter', 1.6, 'beam_divergence', 1.2);
%%

ff = IB.Integrate;
%%
IS.ISopt.mnc = 'HCl_L50_B7_Rw3_rd08_th14.0_D0';
IS.ISopt.beam_diameter = 0.75;
IS.analyze('HR',0,'T',.005,'RD1_margin',2,'r1',3*2.54/2,'r2',3*2.54/2, ...
  'beam_diameter', 0.75, 'beam_divergence', 0);
%%
% This comparison suggests that the ICOS alignment is insensitive to
% beam size and divergence, which is really suprising to me.
load('IB_HCl_L50_B16_Rw3_rd08_th14.0_D12.4_50x100.mat');
ff = IB.Integrate;
close(ff([1 3 4]));
ff1 = ff(2);
figure(ff1);
ax1 = gca;
load('IB_HCl_L50_B7_Rw3_rd08_th14.0_D0.4_50x100.mat');
ff = IB.Integrate;
close(ff([1 3 4]));
ff2 = ff(2);
figure(ff2);
ax2 = gca;
xl = get(ax2,'xlim');
yl = get(ax2,'ylim');
set(ax1,'xlim',xl);
set(ax1,'ylim',yl);
cl = get(ax2,'Clim');
set(ax1,'clim',cl);
%%
IS.res2(4).sel = 1;
IS.res2(7).sel = 0;
for B = 1:5
  for D = 0:3
    IS.ISopt.mnc = sprintf('HCl_L50_B%d_Rw3_D%d',B,D);
    IS.ISopt.beam_diameter = B/10;
    IS.analyze('HR',0,'T',.005,'RD1_margin',2,'r1',3*2.54/2,'r2',3*2.54/2, ...
      'beam_diameter', IS.ISopt.beam_diameter, 'beam_divergence', D);
    IBfile = sprintf('IB_%s.4_50x100.mat', IS.ISopt.mnc);
    load(IBfile);
    ff = IB.Integrate;
    close(ff([1,3,4]));
  end
end
%%
res = collect_results('files','IS_HCl_L50_B*_Rw3_D*.mat');
res = res([res.index]==4);
for i=1:length(res)
  B = sscanf(res(i).mnc, 'HCl_L50_B%d_Rw3_D%d');
  res(i).beam_diameter = B(1)/10;
  res(i).beam_divergence = B(2);
end
%%
scatter([res.beam_diameter],[res.beam_divergence],[],[res.max_pwr]);
xlabel('beam\_diameter');
ylabel('beam\_divergence');
CB = colorbar;
CBL = get(CB,'label');
set(CBL,'string','max\_pwr');
%%
diams = [res.beam_diameter];
divs = [res.beam_divergence];
pwr = [res.max_pwr];
figure;
udiams = unique(diams);
for diam = udiams
  v = [res.beam_diameter] == diam;
  plot(divs(v),pwr(v),'*-');
  hold on;
end
hold off;
xlabel('Divergence degrees');
ylabel('max\_pwr');
legend(num2str(udiams'));
title('Power by beam diameter');
%%
figure;
udivs = unique(divs);
for div = udivs
  v = [res.beam_divergence] == div;
  plot(diams(v),pwr(v),'*-');
  hold on;
end
hold off;
xlabel('Beam Diameter');
ylabel('max\_pwr');
legend(num2str(udivs'));
title('Power by beam divergence');

%%
% Model on-axis alignment with M1 removed
load('IS_HCl_L50_B1_Rw3_D0.mat');
chdir('M2');
B = 1.6;
D = 1.2;
IS.ISopt.mnc = sprintf('HCl_L50_M2_B%d_Rw3_D%.1f',B*10,D);
IS.ISopt.beam_diameter = B;
IS.res2(7).sel = 0;
IS.res2(4).sel = 1;
IS.analyze('HR',0,'T',.005,'RD1_margin',2,'r1',3*2.54/2,'r2',3*2.54/2, ...
  'beam_diameter', IS.ISopt.beam_diameter, 'beam_divergence', D, ...
  'y0', 0, 'z0', 0, 'dy', 0, 'dz', 0,'ICOS_passes',1000);
IBfile = sprintf('IB_%s.4_1000x100.mat', IS.ISopt.mnc);
load(IBfile);
ff = IB.Integrate;
close(ff([1,3,4]));
chdir('..');
%%
% Model on-axis alignment with M1 removed with telescope focused
load('IS_HCl_L50_B1_Rw3_D0.mat');
chdir('M2');
B = 0.8;
D = 0;
IS.ISopt.mnc = sprintf('HCl_L50_M2F_B%d_Rw3_D%.1f',B*10,D);
IS.ISopt.beam_diameter = B;
IS.res2(7).sel = 0;
IS.res2(4).sel = 1;
IS.analyze('HR',0,'T',.005,'RD1_margin',2,'r1',3*2.54/2,'r2',3*2.54/2, ...
  'beam_diameter', IS.ISopt.beam_diameter, 'beam_divergence', D, ...
  'y0', 0, 'z0', 0, 'dy', 0, 'dz', 0,'ICOS_passes',1000);
IBfile = sprintf('IB_%s.4_1000x100.mat', IS.ISopt.mnc);
load(IBfile);
ff = IB.Integrate;
close(ff([1,3,4]));
chdir('..');
%%
% Model on-axis alignment with both mirrors
load('IS_HCl_L50_B1_Rw3_D0.mat');
B = 1.6;
D = 1.2;
IS.ISopt.mnc = sprintf('HCl_L50_OO_B%d_Rw3_D%.1f',B*10,D);
IS.ISopt.beam_diameter = B;
IS.res2(7).sel = 0;
IS.res2(4).sel = 1;
IS.analyze('HR',0,'T',.005,'RD1_margin',2,'r1',3*2.54/2,'r2',3*2.54/2, ...
  'beam_diameter', IS.ISopt.beam_diameter, 'beam_divergence', D, ...
  'y0', 0, 'z0', 0, 'dy', 0, 'dz', 0,'ICOS_passes',50);
IBfile = sprintf('IB_%s.4_50x100.mat', IS.ISopt.mnc);
load(IBfile);
ff = IB.Integrate;
close(ff([1,3,4]));
%%
% Model on-axis alignment with both mirrors
load('IS_HCl_L50_B1_Rw3_D0.mat');
B = 0.8;
D = 0;
IS.ISopt.mnc = sprintf('HCl_L50_OO_B%d_Rw3_D%.1f',B*10,D);
IS.ISopt.beam_diameter = B;
IS.res2(7).sel = 0;
IS.res2(4).sel = 1;
IS.analyze('HR',0,'T',.005,'RD1_margin',2,'r1',3*2.54/2,'r2',3*2.54/2, ...
  'beam_diameter', IS.ISopt.beam_diameter, 'beam_divergence', D, ...
  'y0', 0, 'z0', 0, 'dy', 0, 'dz', 0,'ICOS_passes',50);
IBfile = sprintf('IB_%s.4_50x100.mat', IS.ISopt.mnc);
load(IBfile);
ff = IB.Integrate;
% close(ff([1,3,4]));
