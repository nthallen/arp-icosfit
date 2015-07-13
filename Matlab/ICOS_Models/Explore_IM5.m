%%
P = ICOS_Model5.props;
P.visible = 1;
P.HR = 0;
P.RC = 75;
P.visibility = [0 1 1 1];
P.focus = 1;
P.detector_spacing = 8;
P.ICOS_passes_per_injection = 70;
P.max_rays = 300;
P.y0 = 3.0;
P.injection_scale = .5;
P.dy = 0.097696842105263; % circular
P.dz = 0.057006315789474;
PM = ICOS_Model5(P);

%%
% Check my math regarding propagation of skew and divergence through
% the mirror
[xyz,dxyz,divergence,skew] = PM.M.extract_endpoints_skew(3);
x = 1:length(divergence);
y = ones(size(x));
dcalc = sqrt(sum(xyz(1,[2 3]).^2,2))/P.RC;
scalc = dcalc * sqrt((2*P.RC-P.mirror_spacing)/P.mirror_spacing);
figure; plot(x,divergence,'*',x,dcalc*y); ylabel('divergence');
figure; plot(x,skew,'*',x,scalc*y); ylabel('skew');
%%
[xyz,oxyz] = PM.M.extract_endpoints(4);
dxyz = xyz-oxyz;
dxyz = diag(1./dxyz(:,1)) * dxyz;
yz = oxyz(:,[2 3]); % Using origin, not endpoint
dyz = dxyz(:,[2 3]);
r = sqrt(sum(yz.^2,2));
ryz = diag(1./r) * yz; % unit vector in radial direction
divergence = sum(ryz .* dyz,2);
dcalc = 2.4361 * sqrt(sum(xyz(1,[2 3]).^2,2))/P.RC;
skew = abs(-ryz(:,1).*dyz(:,2)+ryz(:,2).*dyz(:,1));
scalc = sqrt(sum(xyz(1,[2 3]).^2,2))/P.RC * sqrt((2*P.RC-P.mirror_spacing)/P.mirror_spacing);
x = (1:length(divergence))';
y = ones(size(x));
figure; plot(x,divergence./(dcalc*y),'*'); ylabel('divergence'); title('After mirror');
figure; plot(x,skew./(scalc*y),'*'); ylabel('skew'); title('After mirror');

%% Obtain circular alignment
delta=2e-2;
%%
P.visible = false;
P.focus = 0;
P.plot_endpoints = 0;
P.evaluate_endpoints = 3;
PM = ICOS_Model5(P,'dy', P.dy + linspace(-delta,+delta,20),'dz',P.dz + linspace(-delta,delta,20));
PM = PM.clean_results;
figure;
PM.plot_results('eccentricity');
[P.dy,P.dz,eccentricity] = PM.identify_minimum('eccentricity');
delta = delta/5;
%% Graphics to compare eccentricy, max_radius and 
metrics = {'max_radius','eccentricity','overlap'};
for i=1:length(metrics)
  metric = metrics{i};
  figure;
  [C.(metric),H] = contour(PM.Results.(PM.Results.x),PM.Results.(PM.Results.y),PM.Results.(metric),20);
  clabel(C.(metric),H);
  xlabel(PM.Results.x);
  ylabel(PM.Results.y);
  title(strrep(metric,'_','\_'));
  set(gca,'DataAspectRatio',[1 1 1]);
end
%%
C = C.max_radius;
%%
f = figure;
ax = axes('parent',f);
i = 1;
while i < size(C,2)
  %fprintf(1,'Level = %.4f n_pts = %d\n', C(1,i), C(2,i));
  dy = C(1,i+1:i+C(2,i));
  dz = C(2,i+1:i+C(2,i));
  ecc = interp2(PM.Results.(PM.Results.x), PM.Results.(PM.Results.y), PM.Results.eccentricity, dy, dz);
  overlap  = interp2(PM.Results.(PM.Results.x), PM.Results.(PM.Results.y), PM.Results.overlap, dy, dz);
  plot(ax,ecc,overlap);
  hold(ax, 'on');
  i = i + C(2,i) + 1;
end
hold off;
xlabel('eccentricity');
ylabel('overlap');
title('Eccentricity vs Overlap on max\_radius contours');
%%
figure;
scatter(PM.Results.eccentricity(:),PM.Results.overlap(:),5,PM.Results.max_radius(:));
xlabel('eccentricity');
ylabel('overlap');
%%
mr = sort(PM.Results.max_radius(:));
figure; plot(mr);
%% Now move on to focusing
P.focus = 1;
P.visible = 1;
P.Lenses = { 'Lens1' };
P.Lenses = { 'Lens1', 'ZC_PM_25_100' };
P.Lens_Space = [0.2, 6];
P.visibility = [0 0 0];
P.evaluate_endpoints = 4+length(P.Lenses);
PM = ICOS_Model5(P,'detector_spacing', 1:9);
figure; PM.plot_results('max_radius');
[ds,~,max_radius] = PM.identify_minimum('max_radius');
P.detector_spacing = ds;
PM = ICOS_Model5(P);
view(90,0);
figure; PM.M.plot_endpoints(4+length(P.Lenses));
%%
dz0 = P.dz;
%%
P.max_rays = 300;
dzscale = (1:10)/10;
rskew = 0*dzscale;
P.visibility = [];
P.title = 'Scan';
dsrange = P.detector_spacing-2:0.1:P.detector_spacing+2;
ds = zeros(size(dzscale));
maxr = zeros(size(dzscale));
r3 = 2.0;
for i = 1:length(dzscale);
  P.dz = dz0*dzscale(i);
  P.injection_scale = 0.5;
  P.visible = false;
  P.plot_endpoints = 0;
  P.plot_r_angle = 0;
  P.evaluate_endpoints = 3;
  PM = ICOS_Model5(P);
  rmax = PM.Results.max_radius;
  P.injection_scale = P.injection_scale * r3 / rmax;
  P.visible = false;
  P.plot_endpoints = 0;
  P.evaluate_endpoints = 4+length(P.Lenses);
clf;
fill([0 P.D_l/2 P.D_l/2 0 0], [0 0 15 15 0], [0 1 0]);
hold on;
  P.plot_r_angle = 1;
  P.title = sprintf('dzscale: %.1f', dzscale(i));
  PM = ICOS_Model5(P,'detector_spacing',dsrange);
  [P.detector_spacing,~,maxr(i)] = PM.identify_minimum('max_radius');
  if P.detector_spacing == min(dsrange) || P.detector_spacing == max(dsrange)
    error('dsrange does not cover minimum: ds = %.1f, range = [%.1f, %.1f]', ...
      P.detector_spacing, dsrange(1), dsrange(end));
  end
  ds(i) = P.detector_spacing;
  % P.plot_endpoints = 4+length(P.Lenses);
  P.plot_r_angle = 0;
  % PM = ICOS_Model5(P);
%   [xyz,dxyz] = PM.M.extract_endpoints_skew(4+length(P.Lenses));
%   ldxyz = sqrt(sum(dxyz.^2,2));
%   dxyz = diag(1./ldxyz)*dxyz;
%   angle = rad2deg(acos(dxyz(:,1)));
%   r = sqrt(sum(xyz(:,[2 3]).^2,2));
%   plot(r, angle,'*');
%   drawnow; shg;

% This was an earlier investigation
%   [xyz,dxyz,div,skew] = PM.M.extract_endpoints_skew(3);
%   r = sqrt(sum(xyz(:,[2 3]).^2,2));
%   rskew(i) = mean(r.*skew);
%   plot(r,skew,'+r',r,div,'*b');
%   xlim([-1 4]); ylim([-.05 .05]);
%   hold on;
%   drawnow;shg
end
hold off
xlabel('radius at detector');
ylabel('angle of incidence at detector');
ttl = sprintf('%s%s, r3 = %.1f', P.Lenses{1}, ...
  sprintf('+%s', P.Lenses{2:end}), r3);
ttl = strrep(ttl,'+,', ',');
title(strrep(ttl,'_','\_'));
% figure; plot(dzscale,ds,'*');
% xlabel('dzscale');
% ylabel('Detector Spacing');
% figure; plot(dzscale,maxr,'*');
% xlabel('dzscale');
% ylabel('Max Radius');

% xlabel('r');
% ylabel('cm/cm');
% legend('Skew','Divergence','location','southeast');

%% Another approach to looking at skew: abandoning to build another below
P.max_rays = 300;
rs = 2.0;
P.injection_scale = 0.5;
dzscale = (1:5)/10;
dsrange = P.detector_spacing-2:0.1:P.detector_spacing+2;
dzrange = dz0*dzscale;
rskew = 0*dzscale;
P.visibility = [];
P.title = 'Scan';
ds = zeros(size(dzscale));
maxr = zeros(size(dzscale));
P.visible = false;
P.plot_endpoints = 0;
P.plot_r_angle = 0;
P.evaluate_endpoints = 3;
P.detector_spacing = 2.7;
P.dz = dz0*dzscale(1);
PM = ICOS_Model5(P);
P.injection_scale = P.injection_scale * r3 / rmax;
P.evaluate_endpoints = 4+length(P.Lenses);
P.fmt_detector_spacing = 'ds: %.1f';
P.fmt_dz = 'skew: %.3f';
P.plot_endpoints = 4+length(P.Lenses);
P.plot_r_angle = 1;
% PM = ICOS_Model5(P,'dz',dzrange,'detector_spacing',dsrange);
PM = ICOS_Model5(P,'dz',dzrange,'D_dY',linspace(0,.4,10));
% PM = ICOS_Model5(P,'dz',dzrange,'D_l',linspace(0.1,.8,8));
figure; PM.plot_results('total_power');

%% Latest focus investigation
P.Lenses = { 'Lens1', 'ZC_NM_25_100' };
P.Lens_Space = [0.2, 8.5];
P.visible = false;
P.plot_endpoints = 0;
P.plot_r_angle = 0;
P.injection_scale = 0.5;
r3 = 3.0;
dz0 = P.dz;
P.evaluate_endpoints = 3; % mirror 3
fprintf(1,'\n');
for i=1:length(P.Lenses)
  fprintf(1,'%s @ %.1f ', P.Lenses{i}, P.Lens_Space(i));
end
PM = ICOS_Model5(P);
P.injection_scale = P.injection_scale*r3/PM.Results.max_radius;
fprintf(1,'injection_scale = %.1f\n', P.injection_scale);

P.evaluate_endpoints = 4+length(P.Lenses); % detector
PM = ICOS_Model5(P,'detector_spacing', 0.5:.1:1.5); % 8:.2:12); %
PM.plot_results('max_radius');
[P.detector_spacing,~,max_radius] = PM.identify_minimum('max_radius');
fprintf(1,'detector_spacing = %.1f, max_radius at detector: %.2f\n', P.detector_spacing, max_radius);
%%
P.visible = 1;
  P.detector_spacing = 0.1;
  P.evaluate_endpoints = 4+length(P.Lenses);
P.visibility = [0 0 0 0];
PM = ICOS_Model5(P);
P.visible = 0;
%%
dzrange = dz0*(1:5)/10;
tp = zeros(length(dzrange),3); % columns: total power, w/ angle, w/ size
D_dY_range = [0:.05:1.5*max_radius - P.D_l/2 + 0.05];
for i = 1:length(dzrange)
  P.dz = dzrange(i);
  P.D_l = max_radius*2+1;
  P.acceptance_angle = [];
  P.evaluate_endpoints = 4+length(P.Lenses); % detector
  PM = ICOS_Model5(P);
  tp(i,1) = PM.Results.total_power;
  P.acceptance_angle = 15;
  PM = ICOS_Model5(P);
  tp(i,2) = PM.Results.total_power;
  P.D_l = 0.2;
  P.plot_endpoints = 4+length(P.Lenses); % detector
  % PM = ICOS_Model5(P,'D_dY',D_dY_range);
  PM = ICOS_Model5(P,'D_dY',D_dY_range,'detector_spacing', ...
    [.15:.1:.35]);
    %P.detector_spacing+[-.8:.1:0.1]);
  PM.plot_results('total_power'); drawnow; shg;
  P.plot_endpoints = 0;
  [~,~,tp(i,3)] = PM.identify_minimum('-total_power');
  % tp(i,3) = PM.Results.total_power;
  fprintf(1,'%d: TP: %g  +Angle: %.0f %% +Size: %.0f %%\n', ...
    i, tp(i,1), 100*tp(i,2)/tp(i,1), 100*tp(i,3)/tp(i,1));
end
P.dz = dz0;

%%
figure;
plot(dzscale*dz0, rskew,'*');
xlabel('dz');
ylabel('r*skew');

%%
P.visible = 1;
P.visibility = [0 0 0 ];
P.evaluate_endpoints = 7;
P.plot_endpoints = 7;
P.Lenses = {'Lens1', 'ZC_PM_25_25', 'ZC_PM_12_12' };
P.Lens_Space = [0.2, 5, 0.5];
P.detector_spacing = 0.3850;
P.focus = 1;
PM = ICOS_Model5(P);
%%
P.visible = false;
PM = ICOS_Model5(P,'detector_spacing', [0.36:.0025:0.41]);
% PM = PM.clean_results;
PM.plot_results('max_radius');
[ds,~,max_radius] = PM.identify_minimum('max_radius'); % 0.1108
%%
P.detector_spacing = ds;
PM = ICOS_Model5(P);
%%
[xyz,dxyz] = PM.M.extract_endpoints_skew(7);
ldxyz = sqrt(sum(dxyz.^2,2));
dxyz = diag(1./ldxyz)*dxyz;
angle = rad2deg(acos(dxyz(:,1)));
figure; plot(angle,'*');
%%
P.visible = 1;
P.visibility = [0 0 0 ];
P.evaluate_endpoints = 8;
P.plot_endpoints = 8;
P.Lenses = {'Lens1', 'ZC_PM_25_25', 'ZC_PM_12_12', 'ZC_NM_12_12' };
P.Lens_Space = [0.2, 5, 0.5, 0.05];
P.detector_spacing = 0.2460;
P.focus = 1;
PM = ICOS_Model5(P);
%%
P.visible = false;
PM = ICOS_Model5(P,'detector_spacing', [0.23:.0002:0.27]);
% PM = PM.clean_results;
PM.plot_results('max_radius');
[ds,~,max_radius] = PM.identify_minimum('max_radius'); % 0.1195
%%
P.detector_spacing = ds;
PM = ICOS_Model5(P);
%%
[xyz,dxyz] = PM.M.extract_endpoints_skew(7);
ldxyz = sqrt(sum(dxyz.^2,2));
dxyz = diag(1./ldxyz)*dxyz;
angle = rad2deg(acos(dxyz(:,1)));
figure; plot(angle,'*');

%%
P.visible = 1;
P.visibility = [0 0 0 0 0 0 0];
P.evaluate_endpoints = 9;
P.plot_endpoints = 9;
P.LensTypes.ZC_HS_3_3.type = 'negative_meniscus';
P.LensTypes.ZC_HS_3_3.r = 0.2;
P.LensTypes.ZC_HS_3_3.R1 = 0.2;
P.LensTypes.ZC_HS_3_3.R2 = 10.0;
P.LensTypes.ZC_HS_3_3.CT = 0.1;
P.LensTypes.ZC_HS_3_3.EFL = 0.15;
P.Lenses = {'Lens1', 'ZC_PM_25_25', 'ZC_PM_12_12', 'ZC_NM_12_12', 'ZC_HS_3_3' };
P.Lens_Space = [0.2, 5, 0.5, 0.05, 0.098];
P.detector_spacing = 0.2460;
P.focus = 1;
PM = ICOS_Model5(P);
%%
P.visible = false;
PM = ICOS_Model5(P,'detector_spacing', [0.23:.0002:0.27]);
% PM = PM.clean_results;
PM.plot_results('max_radius');
[ds,~,max_radius] = PM.identify_minimum('max_radius'); % 0.1195
%%
P.detector_spacing = ds;
PM = ICOS_Model5(P);
%%
[xyz,dxyz] = PM.M.extract_endpoints_skew(4+length(P.Lenses));
ldxyz = sqrt(sum(dxyz.^2,2));
dxyz = diag(1./ldxyz)*dxyz;
angle = rad2deg(acos(dxyz(:,1)));
r = sqrt(sum(xyz([2 3]).^2,2));
figure; plot(r, angle,'*');

%% Herriot reinjection
% for the 75x75 cylindrical orientation with an radius of 2cm,
% the optimal Herriott mirror (given the choices from Edmund's)
% is a 3" diameter with an 8" FL (Stock #43-582 through 586)
% with a spacing of 12.6 cm
P.visible = 1;
P.visibility = [];
P.focus = 0;
P.injection_scale = 0.7;
P.plot_endpoints = 1;
P.evaluate_endpoints = 3;
P.fmt_HRC = 'Herriott RoC = %.2f';
P.pause = 0;
P.stop_ICOS = 0;
P.HR = 1;
P.HRC = 40.64;
P.herriott_spacing = 12.6;
P.ICOS_passes_per_injection = 63;
P.max_rays = 3000;
%%
P.visible = 1;
P.focus = 1;
PM = ICOS_Model5(P);
%%
PM = ICOS_Model5(P,'HRC',[15.24 30.48 40.64 60.96 91.44 152.4],'herriott_spacing', 7:.2:20);
PM.clean_results;
PM.plot_results('total_power');
%%
P.visible = 0;
P.plot_endpoint = 1;
P.HRC = 40.64;
PM = ICOS_Model5(P,'herriott_spacing',7:.4:20);
%%
PM.plot_results('mean_angle');
%%
PM.plot_results('max_radius');
%%
PM.plot_results('overlap');
%%
PM.plot_results('eccentricity');
%%
PM.plot_results('total_power');
%%
P.HRC = 40.64;
P.herriott_spacing = 12.8;
P.visible = true;
PM = ICOS_Model5(P);

%%
P.HRC = 40.64;
P.herriott_spacing = 10.8;
P.visible = false;
P.plot_endpoints = 1;
P.evaluate_endpoints = 3;
PM = ICOS_Model5(P,'herriott_spacing', [10.8:.2:12.8]);
figure; PM.plot_results('total_power');
%%
P.HRC = 40.64;
P.herriott_spacing = 12.8;
P.visible = false;
P.plot_endpoints = 1;
P.evaluate_endpoints = 3;
PM = ICOS_Model5(P,'injection_scale', [0.5:.1:1]);
figure; PM.plot_results('total_power');

%% Focusing
P.focus = 1;
P.evaluate_endpoints = 5;
PM = ICOS_Model5(P);
%%
P.HR = 0;
P.visibility = [0 0 0];
PM = ICOS_Model5(P);
%%
view(90,0);