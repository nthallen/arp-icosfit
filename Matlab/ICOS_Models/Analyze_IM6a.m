%%
% This configuration is the result of analysis with Explore_IM6a.m
% This design is hypothetical, as we do not currently have the
% second ICOS mirror (ZC-PX-25-500) or the focusing lens
% (ZC-PM-25-38). We also do not have the necessary RIM, although
% the size of this RIM is reasonable.
P = ICOS_Model6.props;
P.visible = 0;
P.HR = 0;
P.visibility = [];
P.focus = 1;
P.y0 = 3.0;
P.z0 = -0.15;
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
P.evaluate_endpoints = -1;
IB = ICOS_beam(@ICOS_Model6, P);
%%
IB.Sample('opt_n', 5, 'ICOS_passes', 10000);
IB.Integrate;
%%
save IB6a_150707_10000x100.mat IB
