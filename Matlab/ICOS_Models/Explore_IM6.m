%%
P = ICOS_Model6.props;
P.visible = 1;
P.HR = 0;
% P.RC = 150; % apparently does not apply
P.visibility = [0 1 1 1];
P.focus = 1;
% P.detector_spacing = 8;
P.ICOS_passes_per_injection = 70;
P.max_rays = 300;
P.R2 = -28.12;
P.y0 = 3.0;
%P.dy = 0.097815789; % circular for -25.4, which isn't available
%P.dz = 0.002078947;
P.dy = 0.097817789473684; % circular for -28.12 which is
P.dz = 0.005696000000000;
P.detector_spacing = 5.105263157894736;
P.Lenses = { 'ZC_PM_25_100' };
P.injection_scale = 2/3;
PM = ICOS_Model6(P);
%%
P.detector_spacing = 2.684210526315789;
P.Lenses = { 'ZC_PM_25_38' };
PM = ICOS_Model6(P);
%%
% Obtain a circular alignment
P.visible = false;
%P.visible = true;
P.focus = 0;
P.plot_endpoints = 3;
P.evaluate_endpoints = 3;
P.injection_scale = .5;
PM = ICOS_Model6(P,'dy', linspace(0.07,0.13,20),'dz',linspace(0,0.03,20));
PM = PM.clean_results;
PM.plot_results('eccentricity');
[P.dy,P.dz,eccentricity] = PM.identify_minimum('eccentricity');
%%
delta=2.5e-3;
%%
P.visible = false;
P.focus = 0;
P.plot_endpoints = 3;
P.evaluate_endpoints = 3;
P.injection_scale = .5;
PM = ICOS_Model6(P,'dy', P.dy + linspace(-delta,+delta,21),'dz',P.dz + linspace(-delta,delta,21));
PM = PM.clean_results;
figure;
PM.plot_results('eccentricity');
[P.dy,P.dz,eccentricity] = PM.identify_minimum('eccentricity');
delta = delta/5;

%% Look at beam spread
P.visible = false;
P.focus = 0;
P.plot_endpoints = 0;
P.evaluate_endpoints = 2;
P.injection_scale = 1;
P.z0 = 0;
%P.injection_scale = 2/3; % 1.0;
%P.z0 = -0.15;
delta = P.beam_diameter/2;
PM = ICOS_Model6(P,'beam_dy', linspace(-delta,+delta,21),'beam_dz',linspace(-delta,delta,21));
%PM = PM.clean_results;
figure;
PM.plot_results('total_power');
%%
figure;
PM.plot_results('overlap');

%% Now adjust focus
P.focus = 1;
P.evaluate_endpoints = 4;
P.visible = 1;
PM = ICOS_Model6(P,'detector_spacing',linspace(2,20,20));
%%
PM.plot_results('max_radius');
%%
[P.detector_spacing,~,max_radius] = PM.identify_minimum('max_radius');
%%
P.Lenses = { 'ZC_PM_25_100' };
PM = ICOS_Model6(P);
%% Now readjust focus
P.focus = 1;
P.evaluate_endpoints = length(P.Lenses) + 4;
P.visible = 0;
PM = ICOS_Model6(P,'detector_spacing',linspace(3,8,20));
PM.plot_results('max_radius');
%%
[P.detector_spacing,~,max_radius] = PM.identify_minimum('max_radius');
%%
P.visible = 1;
P.injection_scale = 1.16;
P.injection_scale = 1.1;
PM = ICOS_Model6(P);
view(50,9);
title(sprintf('R1 = %.0f, R2 = %.1f, L = %.0f', P.R1, P.R2, P.mirror_spacing));
%%
figure; PM.M.plot_endpoints(length(P.Lenses)+4);
ylabel('cm');
title(sprintf('Spot pattern on detector: R1 = %.0f, R2 = %.1f, L = %.0f', P.R1, P.R2, P.mirror_spacing));
%%
figure; PM.M.plot_endpoints(3);
ylabel('cm');
title(sprintf('Spot pattern on mirror 2: R1 = %.0f, R2 = %.1f, L = %.0f', P.R1, P.R2, P.mirror_spacing));
%%
[xyz,dxyz] = PM.M.extract_endpoints_skew(length(P.Lenses)+4);
ldxyz = sqrt(sum(dxyz.^2,2));
dxyz = diag(1./ldxyz)*dxyz;
angle = rad2deg(acos(dxyz(:,1)));
figure; plot(angle,'*');
xlabel('Spot number');
ylabel('degrees');
title(sprintf('Incident angle on detector: R1 = %.0f, R2 = %.1f, L = %.0f', P.R1, P.R2, P.mirror_spacing));
%%
% Enable the Herriott Cell
P.HR = 1;
P.Hr = 8;
P.HRC = 74.84; % 68.5;
P.herriott_spacing = 44.23; % 37.91;
P.focus = 1;
P.evaluate_endpoints = 3;
P.injection_scale = 1; % 1
P.plot_endpoints = 3;
P.visible = 1;
P.visibility = [];
P.ICOS_passes_per_injection = 30;
P.max_rays = 10000;
PM = ICOS_Model6(P);
% PM = ICOS_Model6(P, 'HRC',74:.02:75.5);
%%
P.injection_scale = 2/3;
P.herriott_spacing = 66.34;
P.HRC = 96.939;
P.visible = 1;
P.plot_endpoints = 1;
P.evaluate_endpoints = 1;
PM = ICOS_Model6(P);
%PM = ICOS_Model6(P, 'HRC', 96.93:.001:96.95);
%PM.plot_results('eccentricity');
%%
PM.plot_results('total_power');
%%
PM.plot_results('eccentricity');


view(270,45);