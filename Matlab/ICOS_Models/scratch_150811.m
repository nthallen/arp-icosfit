%% Looking at beam_dy/dz for better power
load('IS_w50_L50.mat')
P = render_model(IS.res2(3)); % R2=-69.8, RL=40.9
%%
P.ICOS_passes_per_injection = 1;
P.visible = 0;
P.focus = 0;
P.evaluate_endpoints = 3;
db = linspace(-0.2, 0.2, 11);
PM = ICOS_Model6(P,'beam_dy',db,'beam_dz',db);
figure;
PM.plot_results('total_power');
% Looks like power for this case would benefit from moving z0=-.05 or a
% bit more. But note that Herriott Cell loss was only 5.6%. ICOS cell
% loss was 28.2% though, so we could pick up more there.
%%
% Rerun analyses