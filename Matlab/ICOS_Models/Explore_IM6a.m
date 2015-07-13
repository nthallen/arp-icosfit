%%
P = ICOS_Model6.props;
P.visible = 1;
P.HR = 0;
P.visibility = [0 1 1 1];
P.focus = 1;
P.ICOS_passes_per_injection = 70;
P.max_rays = 300;
% ZC-PX-25-500 1" dia, RC -69.817 cm
P.r2 = .5*2.54;
P.R2 = -69.817;
P.CT2 = .22;
P.y0 = 3.0;
P.dy = 0.097810000000000;
P.dz = 0.017430000000000;
P.injection_scale = 2/3;
P.detector_spacing = 3.1842;
P.Lenses = { 'ZC_PM_25_38' };
P.Lens_Space = [1.0];
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
P.injection_scale = 2/3;
P.z0 = -.15;
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
P.visible = 0;
P.plot_endpoints = 5;
P.evaluate_endpoints = 5;
%PM = ICOS_Model6(P,'detector_spacing',linspace(2,20,20));
PM = ICOS_Model6(P,'detector_spacing',linspace(6,8,21));
%%
PM.plot_results('max_radius');
%%
[P.detector_spacing,~,max_radius] = PM.identify_minimum('max_radius');
%%
P.Lenses = { 'ZC_PM_25_38' };
P.Lens_Space = [1.0];
PM = ICOS_Model6(P);
%% Now readjust focus
P.focus = 1;
P.evaluate_endpoints = length(P.Lenses) + 4;
P.visible = 0;
PM = ICOS_Model6(P,'detector_spacing',linspace(2.5,3.8,20));
PM.plot_results('max_radius');
%%
[P.detector_spacing,~,max_radius] = PM.identify_minimum('max_radius');
%%
P.visible = 1;
%P.injection_scale = 1.16;
%P.injection_scale = 1.1;
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
P.Hr = 1.5*2.54;
P.HRC = 40.64; % 48.325; % 68.5;
P.herriott_spacing = 18; % 12.6583; % 18; % 37.91;
P.focus = 1;
P.evaluate_endpoints = 1;
P.injection_scale = 2/3;
P.evaluate_endpoints = 1;
P.plot_endpoints = 1;
P.visible = 1;
P.visibility = [1 1 0 0 0];
P.ICOS_passes_per_injection = 30;
P.max_rays = 10000;
PM = ICOS_Model6(P);
%%
P.visible = 0;
P.evaluate_endpoints = 1;
P.plot_endpoints = 1;
PM = ICOS_Model6(P, 'HRC',47.8:.025:48.5);
PM.plot_results('eccentricity');
%%
PM.plot_results('total_power');
%%
PM.plot_results('eccentricity');


view(270,45);