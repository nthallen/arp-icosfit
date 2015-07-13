%% Target practice
P = ICOS_Model5a.props;
P.visible = 1;
P.HR = 1;
P.HRC = 24*2.54;
P.RC = 150;
P.herriott_spacing = 8;
P.visibility = [1 1 1];
P.focus = 0;
P.ICOS_passes_per_injection = 30;
P.max_rays = 1000;
P.y0 = -2.0;
P.injection_scale = 1;
P.dy = -0.031894736842106; % circular
P.dz = 0.023894736842106;
PM = ICOS_Model5a(P);


%% Obtain circular alignment
delta=2e-2;
%%
P.visible = false;
P.focus = 0;
P.plot_endpoints = 0;
P.evaluate_endpoints = 3;
PM = ICOS_Model5a(P,'dy', P.dy + linspace(-delta,+delta,20),'dz',P.dz + linspace(-delta,delta,20));
PM = PM.clean_results;
figure;
PM.plot_results('eccentricity');
[P.dy,P.dz,eccentricity] = PM.identify_minimum('eccentricity');
delta = delta/5;

%% ICOS_Model5a won't play well with this. y0,dy and dz are positions on
% the face of the herriott mirror. To determine optimal herriott mirror
% spacing, we want to use ICOS_Model5, which fixes the position on the
% back of the first ICOS mirror, but we'd need to rotate to get that right.
% P.visible = 0;
% PM = ICOS_Model5a(P,'herriott_spacing',4:.2:8);
% PM = PM.clean_results;
% PM.plot_results('eccentricity');

%%
P.injection_scale = 1.0;
draw_targets(P,1/2.54);
print -dpng -r300 target_2.0.png
%%
P.injection_scale = 1.5;
draw_targets(P,1/2.54);
print -dpng -r300 target_3.0.png
