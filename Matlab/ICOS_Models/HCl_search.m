%%
% HCL_search
rd = 0.08; % Define the target radius at the detector
th = 13; % Define the target angle from normal at the detector
L = [50 100]; % Cavity length
R2_lim = [-1000 0];

R1 = [50 75 100 200 300 400 500 750 1000];
IM2 = ispcatalog;
IR2 = unique([IM2.R_cm]);
if ~isempty(R2_lim)
  IR2 = IR2(IR2 > R2_lim(1) & IR2 < R2_lim(2));
end
R2 = [IR2 R1];
HM1 = ed_rim_catalog;
RR1 = unique([HM1.R_cm]);
%%
B = 0.1;
mnc = sprintf('HCl_B%d_rd%02d_th%.1f', round(B*10), floor(rd*100), th);
SR = ICOS_sr_search('mnc',mnc,'L',L,'B',B,'Rw1',0.2,'R1',R1,'R2',R2,...
  'RR1',RR1,'rd',rd,'th',th, 'ICOS_margin', 0.2, 'RD1_margin', 2.0);
SR.enumerate;
%%
SR.design_tolerance;
SR.build_tolerance;
%%
ddBL = diff([SR.Summary.dBL])';
r_max = [SR.Summary.r_max]';
mirror_resolution = 2.54/2;
mirror_margin = SR.SRopt.B * 0.75;
mirror_lip = 2.54/8; % 1/8"
r_mirror = ceil((r_max+mirror_margin+mirror_lip)/mirror_resolution) * ...
  mirror_resolution - mirror_lip;
% r_mirror = r_max+mirror_margin+mirror_lip; % This is theoretical min
L = [SR.Summary.L]';
Volume = pi*L.*(r_mirror-mirror_lip).^2*1e-3; % liters
%%
% Now pick configurations based on ddBL, V
ddBLok = ddBL > 0.2;
N_configs = 20;
[~,VI] = sort(Volume);
Vlow = sum(cumsum(ddBLok(VI))<=N_configs);
SR.select(VI(1:Vlow));
SR.deselect(find(~ddBLok));
% Skip ahead to savefile
%%
for i=1:5
  SR.explore_build_tolerance(i);
end
%%
v = find([SR.Summary.RLmin] > 28);
SR.deselect(v);
%%
dBL = [SR.Summary.dBL];
ddBL = diff(dBL);
v = find(ddBL < 0.2);
SR.deselect(v);
%%
r2r1 = [SR.Summary.r2]./[SR.Summary.r1];
sel = find([SR.Summary.sel]);
fprintf(1,'Minimum r2/r1 selected = %.2f\n', min(r2r1(sel)));
%%
sel = [SR.Summary.sel];
v = find(~sel & (r2r1 < 0.4) & (ddBL > 0.2));
fprintf(1,'Found %d new configurations\n', length(v));
%%
RLmin = [SR.Summary.RLmin];
L = [SR.Summary.L];
v = find(~sel & (r2r1 < 0.4) & (ddBL > 0.2) & (RLmin < 40));
fprintf(1,'Found %d new configurations\n', length(v));
%%
v = find(~sel & (r2r1 < 0.4) & (ddBL > 0.2) & (RLmin < 40) & (L<65.75));
fprintf(1,'Found %d new configurations\n', length(v));
%%
SR.select(v);
%%
sel = [SR.Summary.sel]>0;
RLmin = [SR.Summary.RLmin];
RR1 = [SR.Summary.RR1];
L = [SR.Summary.L];
r2r1 = [SR.Summary.r2]./[SR.Summary.r1];
%%
figure;
plot(L,RLmin,'.',L(sel),RLmin(sel),'or');
xlabel('L');
ylabel('RLmin');
%%
figure;
plot(r2r1,L+RLmin,'.',r2r1(sel),L(sel)+RLmin(sel),'or');
xlabel('r_2/r_1');
ylabel('L+RL');

%%
SR.savefile;
%%
SR.focus('focus_visible',0);
%%
results = dir('IS_HCl*.mat');
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
pattern = 'IS_HCl_rd08_th13.0*.mat';
%%
res = collect_results('files',pattern,'R1',[],'D1',[],'sel',1,'exclude','NH','exclude','max_pwr');
%%
% Discard the IS files with no configurations
v = find(configs==0);
for vi = v'
  fprintf(1,'Deleting %s\n', results(vi).name);
  delete(results(vi).name);
end
fprintf(1,'Remember to recalculate results and focuses\n');
%%
% Whoa! I only got successful focus on 10 out of 75 configurations
% What was the issue?
s2 = zeros(length(results),1);
r2 = zeros(length(results),1);
d2 = zeros(length(results),1);
d3 = zeros(length(results),1);
for i=1:length(results)
  load(results(i).name);
  s2(i) = IS.res1.s2;
  r2(i) = IS.res1.r2;
  d2(i) = IS.res1.d2;
  d3(i) = d2(i)*IS.res1.n;
end
th = atand(sqrt(s2.^2 + d3.^2));
v = focuses>0;
figure;
%%
rdtanth = .08 * tand(13);
plot(th(v),s2(v).*r2(v),'*b',th(~v),s2(~v).*r2(~v),'*r',minmax(th'),rdtanth*[1 1],'g');
xlabel('theta'); ylabel('s_2r_2'); shg;
%%
plot(s2(v),r2(v),'*b',s2(~v),r2(~v),'*r');
xlabel('s2'); ylabel('r2');
shg;
%%
plot(s2(v),d3(v),'*b',s2(~v),d3(~v),'*r');
xlabel('s2'); ylabel('d3');
shg;
%%
plot(d3(v),s2(v).*r2(v),'*b',d3(~v),s2(~v).*r2(~v),'*r',minmax(d3'),rdtanth*[1 1],'g');
xlabel('d_3'); ylabel('s_2r_2'); shg;
%%
results = dir('IS_HCl*.mat');
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
% Manually review focus
results = dir('IS_HCl*.mat');
for i=1:length(results)
  if focuses(i) > 0
    load(results(i).name);
    IS.explore_focus;
    IS.savefile;
  end
end
%%
% Automatically choose focuses
results = dir('IS_HCl*.mat');
Max_Foci = 2;
for i=1:length(results)
  if focuses(i) > 0
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
end
%%
res = collect_results('files','IS_HCl_rd08_th13.0*.mat','sel',1,'D1',[],'exclude','NH','exclude','max_pwr');
for i=1:length(res)
  radius = res(i).D1*2.54/2;
  res(i).Volume = res(i).L*pi*radius^2*1e-3;
end
%%
results = dir('IS_HCl*.mat');
for i=1:length(results)
  if focuses(i) > 0
    load(results(i).name);
    mnc = IS.ISopt.mnc;
    SRi = strrep(mnc,[SR.SRopt.mnc '.'],'');
    len = min(strfind(SRi,'_')) - 1;
    SRi = str2num(SRi(1:len));
    r_max = SR.Summary(SRi).r_max;
    r1 = SR.Summary(SRi).r1;
    injection_scale = r_max/r1;
    fprintf(1,'%d: injection_scale = %f\n', i, injection_scale);
    IS.analyze('HR',0,'injection_scale',injection_scale); % Analysis w/o Herriott cell
    close all;
  end
end
%%
res = collect_results('files','IS_HCl_rd08_th13.0*.mat');
%%
max_pwr = [res.max_pwr];
r1 = [res.r1];
L = [res.L];
mirror_resolution = 2.54/2;
mirror_margin = SR.SRopt.B * 0.75;
mirror_lip = 2.54/8; % 1/8"
r_mirror = ceil((r1+mirror_margin+mirror_lip)/mirror_resolution) * ...
  mirror_resolution;
% r_mirror = r_max+mirror_margin+mirror_lip; % This is theoretical min
Volume = pi*L.*(r_mirror-mirror_lip).^2*1e-3; % liters
overlap = [res.overlap];
Rw1 = [res.Rw1];
%%
plot(max_pwr,overlap,'.');
%%
plot(max_pwr,Volume,'.');
%%
plot(max_pwr,Rw1,'.');
%%
results = dir('IS_HCl_*.mat');
for i=1:length(results)
  if focuses(i) > 0
    load(results(i).name);
    mnc = strrep(IS.ISopt.mnc,'HCl_','HCl.H_');
    IS.ISopt.mnc = mnc;
    IS.savefile;
    IS.analyze('RD1_margin',2);
    close all;
    % SRi, R1, R2, RR1, m, k
    %A = sscanf(IS.ISopt.mnc,[SR.SRopt.mnc '.%d_%d_%d_%d_%d.%d']);
  end
end
%%
res = collect_results('files','IS_HCl.H_rd08_th13.0*.mat','Rr1',[], ...
  'r2',[],'D1',[],'Lens_Space',[],'detector_spacing',[]);
%%
for i=1:length(res)
  mnc = res(i).mnc;
  % SR#, R1, R2, RR1, m, k
  A = sscanf(mnc,'HCl.H_rd08_th13.0.%d_%d_%d_%d_%d.%d');
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
% eNH does not seem correlated with Rw1. Where is the power going?
% quick_plot(res,'Rw1','eNH'); quick_plot(res,'Rw1');
% shows that Rw1 = 0.4137, and this is for samples 19, 20. 20 has the
% lower eNH
resno = 11;
IBfile = sprintf('IB_%s.%d_50x100.mat', res(resno).mnc, res(resno).index);
load(IBfile);
mnc = strrep(res(resno).mnc,'HCl.H','HCl.H1');
P = IB.P;
P.stop_ICOS = 1;
P.focus = 0;
% P.beam_diameter = 0.2;
%% Animate optic 1
IB = ICOS_beam(@ICOS_Model6,P);
IB.Sample('mnc', mnc, 'opt_n',1,'n_optics',3,'ICOS_passes',1);
IB.Animate('IPass',0);
%% Animate optic 2
IB = ICOS_beam(@ICOS_Model6,P);
IB.Sample('mnc', mnc, 'opt_n',2,'n_optics',3,'ICOS_passes',1);
IB.Animate('IPass',0);
