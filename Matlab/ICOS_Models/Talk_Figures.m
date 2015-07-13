%%
P = ICOS_Model5.props;
P.visible = 1;
P.HR = 0;
P.herriott_spacing = 2;
P.RC = 75;
P.visibility = [0 1 1 1];
P.focus = 1;
P.detector_spacing = 8;
P.ICOS_passes_per_injection = 70;
P.max_rays = 300;
P.y0 = 3.0;
P.injection_scale = 1;
P.dy = 0.097696842105263; % circular
P.dz = 0.057006315789474;
PM = ICOS_Model5(P);
%%
P.evaluate_endpoints = 5;
P.visibility = [0 0 0 ];
PM = ICOS_Model5(P,'detector_spacing',9:.1:11);
[P.detector_spacing,~] = PM.identify_minimum('max_radius');
%%
PM.plot_results('max_radius');
%%
P.visibility = [0 1];
PM = ICOS_Model5(P);
%%
P.visibility = [0 0 0 0];
PM = ICOS_Model5(P);
%%
view(90,0);
%%
P.visibility = [0 0 0];
P.evaluate_endpoints = 0;
P.plot_endpoints = 0;
P.focus = 1;
PM = ICOS_Model5(P);

%% Can going to a more eccentric alignment help?
P.max_rays = 2000;
P.ICOS_passes_per_injection = 2000;
P.focus = 0;
P.visible = 0;
P.plot_endpoints = 3;
P.evaluate_endpoints = 3;
dzf = 1:-.1:.01;
PM = ICOS_Model5(P,'dz', P.dz * dzf);
