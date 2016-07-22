%%
% HCL_search_B
% This is very much like HCl_search, but without most of the exploration
% code. It just runs straight through.
%
% The current adaptation is to look at whether there are shorter configs
% that we were excluding arbitrarily.
% mnc starts with 'HCl_B4_rd...'
rd = 0.08; % Define the target radius at the detector
th = 14; % Define the target angle from normal at the detector
L = 50; % Cavity length
% R2_lim = [-1000 0];
R2_lim = [];

% R1 = [25 50 75 100 200 300 400 500 750 1000];
R1 = 75;
IM2 = ispcatalog;
IR2 = unique([IM2.R_cm]);
if ~isempty(R2_lim)
  IR2 = IR2(IR2 > R2_lim(1) & IR2 < R2_lim(2));
end
% R2 = [IR2 R1];
R2 = R1;
HM1 = ed_rim_catalog;
RR1 = unique([HM1.R_cm]);
%
B = 0.4;
ICOS_margin = 0.2;
mnc = sprintf('HCl_L50_B%d_rd%02d_th%.1f', round(B*10), floor(rd*100), th);
%%
IS = ICOS_search('mnc',mnc,'L',L,'R1',R1,'R2',R2,'Rw1',B/2,'RL_lim',[5, 50]);
IS.search_ICOS_RIM;
%%
IS.search_focus2('focus_visible', 0, 'max_lenses', 2);
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
IS.analyze('RD1_margin',2);

%%
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
