%%
% Telescope diagnostics
%
% This first section starts with zero divergence and just demonstrates
% that as originally specified, the output beam is somewhat convergent.
P = Telescope.props;
P.visible = 0;
P.evaluate_endpoints = 3;
P.beam_n_r = 1;
P.beam_n_th = 12;
P.L1_Space = 17; % As measured 8/15/16
P.L2_Space = 6.237; % As measured 8/15/16
% P.D_Space = 37; % As measured 8/15/16
PM = Telescope(P,'beam_diameter',[0.1 0.2 0.3 0.4 0.5]);
PM.plot_results('max_divergence');
%%
% Just visualize the last optic
P.visible = 1;
P.beam_diameter = 0.1;
P.visibility = [0 0];
PM = Telescope(P);
%%
% This shows that for parallel inputs, divergence stays pretty small
% at different L2_Space values. Zero is near 6.096 (for 0.4 beam at 1 cm)
P.visible = 0;
P.evaluate_endpoints = 3;
P.beam_n_r = 1;
PM = Telescope(P,'L2_Space',[6.085:.001:6.15]);
PM.plot_results('max_divergence');
grid;
%%
% Now let's look at how this changes when divergence is added
% Output divergence increases at 4*input divergence.
P = Telescope.props;
P.visible = 0;
P.evaluate_endpoints = 3;
P.beam_n_r = 1;
P.beam_n_th = 12;
P.L2_Space = 6.096;
PM = Telescope(P,'beam_divergence',[0:.02:0.3]);
PM.plot_results('max_divergence');
hold on;
plot(PM.Results.beam_divergence,PM.Results.beam_divergence*4);
grid;
legend('output divergence','input divergence * 4','location','southeast');
% Confirms that divergence scales up by 4
%%
PM.plot_results('max_radius');
% this shows beam radius is actually 1/4 and improves slighly with divergence
% That makes sense because in the input divergence pushes the focal point
% out, making the radius smaller at L2, but that means the angular
% correction is less.
%%
% Now let's see what happens if we try to improve the focus:
% It appears we should be able to reduce the divergence by
% pushing L2 out to 6.8 (from ~6.1) without significantly increasing
% the output beam diameter
P = Telescope.props;
P.visible = 0;
P.evaluate_endpoints = 3;
P.L2_Space = 6.096;
P.beam_diameter = 0.4;
P.beam_divergence = 0.3; % degrees
P.beam_n_r = 1;
P.beam_n_th = 12;
dL2S = linspace(0,1,50);
PM = Telescope(P,'L2_Space', P.L2_Space + dL2S);
figure;
PM.plot_results('max_divergence');
grid;
figure;
PM.plot_results('max_radius');
grid;
%%
close all
%%
% This is just a test to verify that we aren't diverging
% Just push the detector out and verify the max_radius does not increase
% appreciably.
P = Telescope.props;
P.visible = 0;
P.evaluate_endpoints = 3;
P.beam_diameter = 0.4;
P.beam_divergence = 0.3; % degrees
P.beam_n_r = 1;
P.beam_n_th = 12;
P.L2_Space = 6.8; % Matched to 0.3
PM = Telescope(P,'D_Space',1:20);
PM.plot_results('max_radius');
%%
PM.plot_results('max_divergence');
%%
% Now let's try to model the HCl system as currently configured, with
% 18 cm between the laser and L1, L2_Space of 6.35 cm and ~20 cm before M1
P = Telescope.props;
P.visible = 1;
P.evaluate_endpoints = 3;
P.beam_diameter = 0.4;
P.beam_divergence = 0.3; % degrees
P.beam_n_r = 1;
P.beam_n_th = 12;
P.L1_Space = 17;
P.L2_Space = 6.24; % As measured 8/15/16
P.D_Space = 37; % Actually it's more like 37, but hard to visualize
PM = Telescope(P);
view(0,0);
%%
% Now let's find an optimal L2_Space
P.beam_n_r = 1;
P.visible = 0;
dL2S = linspace(0,1,50);
PM = Telescope(P,'L2_Space', P.L2_Space + dL2S);
figure;
PM.plot_results('max_divergence');
grid;
figure;
PM.plot_results('max_radius');
grid;
%%
% I get 6.78: Now let's see how this looks:
P.L2_Space = 6.78;
P.visible = 1;
PM = Telescope(P);
view(0,0);

%%
% Test ICOS_beam verison of Telescope using the current L2_Space:
P = Telescope.props;
P.visible = 0;
P.evaluate_endpoints = 3;
P.beam_diameter = 0.4;
P.L1_Space = 18;
P.L2_Space = 6.24;
P.D_Space = 37;
IB = ICOS_beam(@Telescope, P);
IB.Sample('beam_samples', 500, ...
            'opt_n', 3, ...
            'n_optics', 3, 'Track_Power', 0, ...
            'beam_divergence', 0.3, ...
            'mnc', 'Tele');
%
ff = IB.Integrate;
close(ff([1,3,4]));
ff1 = ff(2);
%%
% And the optimized L2_Space:
P.L2_Space = 6.78;
IB = ICOS_beam(@Telescope, P);
IB.Sample('beam_samples', 500, ...
            'opt_n', 3, ...
            'n_optics', 3, 'Track_Power', 0, ...
            'beam_divergence', 0.3, ...
            'mnc', 'Tele');
ff = IB.Integrate;
close(ff([1,3,4]));
ff2 = ff(2);
figure(ff2);
cl = get(gca,'Clim');
figure(ff1);
set(gca,'Clim',cl);
