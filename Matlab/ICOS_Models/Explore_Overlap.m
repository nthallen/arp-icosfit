%% Overlap
P = ICOS_Model5.props;
P.visible = 0;
P.HR = 0;
P.RC = 75;
P.visibility = [0 1 1 1];
P.focus = 0;
P.detector_spacing = 8;
P.ICOS_passes_per_injection = 70;
P.max_rays = 300;
P.y0 = 3.0;
P.injection_scale = .5;
P.dy = 0.097696842105263; % circular
P.dz = 0.057006315789474;
%%
MS = 52:.02:55;
overlap = zeros(size(MS));
eccentricity = zeros(size(MS));
max_radius = zeros(size(MS));
delta=2e-2;
for i=1:length(MS)
  P.mirror_spacing = MS(i);
  while delta > 3e-3
    P.visible = false;
    P.focus = 0;
    P.plot_endpoints = 0;
    P.evaluate_endpoints = 3;
    PM = ICOS_Model5(P,'dy', P.dy + linspace(-delta,+delta,20),'dz',P.dz + linspace(-delta,delta,20));
    PM = PM.clean_results;
    % figure;
    % PM.plot_results('eccentricity');
    fprintf(1,'i=%d, MS=%f, delta=%g\n', i, P.mirror_spacing, delta);
    drawnow;
    [P.dy,P.dz,eccentricity] = PM.identify_minimum('eccentricity');
    delta = delta/5;
  end
  delta = delta*5;
  PM = ICOS_Model5(P);
  overlap(i) = PM.Results.overlap;
  eccentricity(i) = PM.Results.eccentricity;
  max_radius(i) = PM.Results.max_radius;
end
figure; plot(MS,overlap);
xlabel('Mirror Spacing');
ylabel('Overlap');
figure; plot(MS,eccentricity);
xlabel('Mirror Spacing');
ylabel('Eccentricity');
figure; plot(MS,max_radius);
xlabel('Mirror Spacing');
ylabel('max radius');
%%
%save overlap2.mat MS overlap eccentricity max_radius
