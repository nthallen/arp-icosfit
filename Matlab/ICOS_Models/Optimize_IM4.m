%%
cd C:\Data\ES96\ICOS\3D
P = ICOS_Model4.props;

% Formatting
P.visible = false;
% P.avifile = 'L2_sweep.avi';
P.plot_endpoints = 0;
P.pause = 0;
P.title = 'Scan';
P.fmt_L2_X = 'L2\\_X = %.1f cm';
P.fmt_L1_R2 = 'L1\\_R2 = %.2f cm';

P.herriott_spacing = 10; % Before the first ICOS mirror
P.HRC = 6*2.54; % Herriott radius of curvature
P.HR = 0;
P.injection_scale = 1;
P.RC = 150;
P.y0 = P.Hr/2; % Location of Herriott hole
P.dy = .0344; % optimal for some other parameters
P.dz = 0.0328; % optimal for some other parameters
P.max_rays = 200;
P.ICOS_passes_per_injection = 70;
%%
status_fig = figure;
n_L = 20;
Lvals = linspace(50,46.5,n_L);
z = zeros(n_L,1);
Lscan = struct('mirror_spacing',z,'dy',z,'dz',z,'overlap',z,'eccentricity',z);
for i=1:n_L;
  P.visible = false;
  P.focus = 0;
  P.evaluate_endpoints = 3;
  P.mirror_spacing = Lvals(i);
  PM = ICOS_Model4(P,'dy', linspace(0,0.07,20),'dz',linspace(0,0.05,20));
  PM = PM.clean_results;
  P.visible = false;
  [dy,dz,~,ri,ci] = PM.identify_minimum('eccentricity');
  % We can get more precise by doing a second pass:
  PM = ICOS_Model4(P,'dy', dy + linspace(-3.5e-3,+3.5e-3,20),'dz',dz + linspace(-2.5e-3,2.5e-3,20));
  PM = PM.clean_results;
  PM.plot_results('eccentricity');
  [dy,dz,eccentricity] = PM.identify_minimum('eccentricity');
  PM.plot_results('eccentricity');
  title(sprintf('L = %.2f', P.mirror_spacing));
  drawnow;
  % PM.plot_results('overlap');

  % Now save mirror_spacing, dy, dz, overlap, eccentricity
  Lscan.mirror_spacing(i) = P.mirror_spacing;
  Lscan.dy(i) = dy;
  Lscan.dz(i) = dz;
  Lscan.eccentricity(i) = eccentricity;
  Lscan.overlap(i) = PM.Results.overlap(ri,ci);
  fprintf(1,'L: %.2f overlap: %.1f eccentricity: %.2f\n', ...
    P.mirror_spacing, Lscan.overlap(i), eccentricity);
end
figure;
plot(Lscan.mirror_spacing, Lscan.overlap);
%%
% Now pick the mirror_spacing with the minimum overlap, optimal angle
overlap = min(Lscan.overlap);
min_i = find(Lscan.overlap == overlap,1);
P.mirror_spacing = Lscan.mirror_spacing(min_i);
P.dy = Lscan.dy(min_i);
P.dz = Lscan.dz(min_i);
figure(status_fig);
%%
% And switch to the Herriott alignment
P.visible = false;
P.plot_endpoints = 3;
P.evaluate_endpoints = 1;
P.title = 'Scan';
P.fmt_herriott_spacing = 'RIM\\_X = %.2f cm';
P.fmt_HRC = 'RIM\\_RC = %.2f cm';
P.HR = 1;
P.stop_ICOS = true;
P.max_rays = 200;
HRC_vals = [3 6 8 12 18 30] * 2.54 * 2; % focal length in inches
% HRC_vals = [15.24 30.48 40.64 60.96 91.44 152.4];
PM = ICOS_Model4(P,'HRC', HRC_vals, 'herriott_spacing',3:.2:15);
%%
figure;
PM = PM.clean_results;
PM.plot_results('eccentricity');
title(sprintf('RC: %.0f L: %.0f', P.RC, P.mirror_spacing));
figure;
PM.plot_results('n_rays');
figure;
PM.plot_results('total_power');
%%
[HRC, herriott_spacing,eccentricity,ri,ci] = PM.identify_minimum('-total_power');
HFL = HRC / (2.54*2);
fprintf(1,'Optimal values:\n  RC: %.2f cm\n  L: %.2f cm\n  HFL: %.2f in\n  HL: %.2f cm\n  y0: %.3f\n  dy: %.4f\n  dz: %.4f\n', ...
  P.RC, P.mirror_spacing, HFL, herriott_spacing, P.y0, P.dy, P.dz);
%%
P.HRC = HRC;
P.herriott_spacing = herriott_spacing;
P.stop_ICOS = false;
P.visible = true;
P.max_rays = 1250;
P.ICOS_passes_per_injection = 40;
PM = ICOS_Model4(P);
%%
% Now save the results for later analysis
save Lscan150.mat Lscan PM P
%%
figure;
PM.M.plot_endpoints(1);
%%
figure;
PM.M.plot_endpoints(3);

