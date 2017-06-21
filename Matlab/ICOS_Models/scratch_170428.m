%%
% White cell exercises
cd C:\Users\nort.ARP\Documents\SW\arp-icosfit\Matlab\ICOS_Models
%%
P = WhiteCell.props(25,5);
P.visible = 1;
P.beam_samples = 100;
% P.max_rays = P.beam_samples*8*P.N4;
P.beam_divergence = 6.2; % degrees
%P.M0_Ddz = .2/25; % tilt the primary mirror
% P.M0_roc = 27;
P.Cell_Length = 25;

P.Lens_dx = 3.0; % Distance from the housing inner wall
P.LED_dx = 0.3; % Distance from the lens

% % 6mm x 9mm FL
% P.Lens_r = 0.5;
% P.Lens_CT = 0.27;
% P.Lens_ROC = 0.78;
% % 12.5mm x 50mm FL
% P.Lens_r = 1.25;
% P.Lens_CT = .566;
% P.Lens_ROC = 4.494;

% Leaving the rez here
P.Lens_r = 0.5;
P.Lens_type = 'double_convex';
P.Lens_ROC = 3.0;

% Results from optimization:
P.Lens_ROC = 1.9;

% Attempt to use Edmunds part 12mm x 18mm FL plano convex
P.Lens_r = 0.6;
P.Lens_ROC = 0.825;
P.Lens_CT = 0.400;
P.Lens_type = 'plano_convex';
% 12mm x 20mm double convex
P.Lens_r = 0.6;
P.Lens_ROC = 1.763;
P.Lens_CT = 0.437;
P.Lens_EFL = 2.0;
P.Lens_type = 'double_convex';
% 12mm x 24mm double convex #49-253
P.Lens_r = 0.6;
P.Lens_ROC = 2.151;
P.Lens_CT = 0.310;
P.Lens_EFL = 2.4;
P.Lens_type = 'double_convex';

% P.visibility = [1 1 0 0];
P.evaluate_endpoints = 0;
%%
P.rng_state = rng;
%%
PM = WhiteCell(P);
%xlim([-2 2]);
%%
P.beam_samples = 1000;
P.visible = 0;
P.rng_state = rng;
PM = WhiteCell(P);
%%
% Look at power per pass
opt_n = PM.evaluate_passes;
max_pass = length(opt_n);
power = zeros(max_pass,1);
spot_dia = power;
for i=1:max_pass
  P.evaluate_endpoints = i;
  Res = PM.evaluate_endpoints(P);
  power(i) = Res.Ipower;
  spot_dia(i) = Res.spot_dia;
end
loss = 1 - power./[1; power(1:end-1)];
%
figure;
plot(power,'*');
xlabel('Pass');
ylabel('Power');
title(sprintf('M0\\_Ddz = %.4f divergence = %.1f', P.M0_Ddz, P.beam_divergence));
%%
figure;
passes = 1:max_pass;
plot(passes,spot_dia,'*r',[1,max_pass],[1,1]*0.25*2.54/2,[1,max_pass],[1,1]*P.M1_r);
xlabel('Pass');
ylabel('Spot Diameter cm');
title(sprintf('M0\\_Ddz = %.4f divergence = %.1f', P.M0_Ddz, P.beam_divergence));
%%
figure;
plot(loss*100,'*r');
xlabel('Pass');
ylabel('% Loss');
title(sprintf('M0\\_Ddz = %.4f divergence = %.1f', P.M0_Ddz, P.beam_divergence));
%% and draw the figure with individual passes
i = 1;
%%
PM.M.redraw(find(PM.M.User.pass==i)); %#ok<FNDSB>
title(sprintf('Pass %d',i));
shg;
if i == max(PM.M.User.pass)
  i = 1;
else
  i = i+1;
end
%%
P.visible = 0;
P.evaluate_endpoints = 3;
PM = WhiteCell(P,'Lens_dx', [0.1:0.1:1.0]);
%%
P.evaluate_endpoints = 6;
PM6 = WhiteCell(P,'Lens_dx', [0.1:0.1:1.0]);
%%
P.evaluate_endpoints = 1;
PM1 = WhiteCell(P,'Lens_dx', [0.1:0.1:1.0]);
%%
P.pause=0;
P.view = [-90,0];
PM = WhiteCell(P,'M2_fy', -(0:.01:0.5));
%%
% This is obsolete, as I've abandoned ICOS_beam for this analysis
P.visible = 0;
passes = [];
sizes = [];
figure;
for opt_n = 2:4
  IB = ICOS_beam(@WhiteCell,P,'n_optics',4,'beam_samples',100,'beam_divergence',1);
  IB.Sample('opt_n', opt_n);
  passvals = unique(IB.Res.NPass(1:IB.Res.N,1));
  psizes = zeros(length(passvals),1);
  for i = 1:length(passvals)
    V = IB.Res.NPass(1:IB.Res.N,1)==passvals(i);
    psizes(i) = 2*sqrt(var(IB.Res.E(V,2)) + var(IB.Res.E(V,3)));
  end
  if opt_n == 2
    passvals = passvals * 2 - 2;
  else
    passvals = passvals * 2 - 1;
  end
  plot(passvals,psizes,'*');
  hold on;
  passes = [passes passvals'];
  sizes(passvals) = psizes;
end
hold off;
legend('M0','F1','F2');
ylabel('Spot Diameter cm');
xlabel('Pass Number');
passes = unique(passes);
%%
figure;
plot(passes, sizes(passes), '*');
ylabel('Spot Diameter cm');
xlabel('Pass Number');
title(sprintf('Spots on optic %d', IB.IBP.opt_n));
%%
% Try using evaluate_endpoints
P = WhiteCell.props(25,2);
P.visible = 0;
P.beam_samples = 100;
P.beam_divergence = 1; % degree
% P.M0_roc = 27;
P.Cell_Length = 25;
P.evaluate_endpoints = 8;
P.rng_state = rng;
% P.visibility = [1 1 0 0];
%%
PM = WhiteCell(P,'Cell_Length',23:.1:27,'beam_divergence',1:.1:2);
%%
PM = WhiteCell(P,'beam_divergence',0:.1:2);
%%
PM = WhiteCell(P,'M0_Ddy',(0:.1:1)*.2/25,'M0_Ddz',(0:.1:1)*.2/25);
%%
PM = WhiteCell(P,'M0_Ddz',(0:.1:1)*.2/25,'beam_divergence',0:.2:2);
%%
P.beam_divergence = 1;
P.M0_Ddz = .2/25;
PM = WhiteCell(P,'M1_Ddz',(-1:.04:1)*0.5);
%%
figure;
PM.plot_results('Ipower');
%%
figure;
PM.plot_results('Opower');
%%
figure;
PM.plot_results('spot_dia');
%%
%%
% Try using evaluate_endpoints to investigate M0 tilt
P = WhiteCell.props(25,2);
P.visible = 0;
P.beam_samples = 100;
P.beam_divergence = 2; % degree
P.M0_Ddz = .2/25; % tilt the primary mirror
%P.dzdz = .05;
% P.M0_roc = 27;
P.Cell_Length = 25;
P.evaluate_endpoints = 8;
P.rng_state = rng;
% See if tilting M1 can correct for the power loss
P.evaluate_endpoints = 8;
%%
P.visible = 1;
PM = WhiteCell(P);
%%
PM = WhiteCell(P,'M2_Ddz',.1 * (-1:.1:1));
%%
PM = WhiteCell(P,'M1_dy',P.M1_dy + .1 * (-1:.1:1));
%%
PM = WhiteCell(P,'dzdz', .03 * (-1:.05:1),'M1_dy',P.M1_dy + .1 * (-1:.1:1));
%%
figure;
PM.plot_results('Ipower');
%%
figure;
PM.plot_results('Opower');
%%
figure;
PM.plot_results('spot_dia');
%%
% Play with opt_model_p and evaluate endpoints
% P.evaluate_endpoits specifies which *pass* to look at, rather than
% which optic.
%  Try moving the LED forward. The minimum x position at M0 is -P.M0_CT
% The lens has P.Lens_CT, and we have P.Lens_dx before the LED, which is at
% LED_x. So the maximum value for LED_x is -P.M0_CT-P.Lens_CT-P.Lens_dx
LED_x_max = -P.M0_CT-P.Lens_CT-P.Lens_dx;
LED_x_max = -2;
LED_x_min = -4;
P.visible = 0;
P.beam_samples = 100;
P.evaluate_endpoints = 23; % 3 is M1, 23 is detector
P.rng_state = rng;
PM = WhiteCell(P, 'LED_x', LED_x_min + [0:.1:1]*(LED_x_max-LED_x_min), ...
  'Lens_ROC', [1.5:0.1:2.5]);
%%
figure;
PM.plot_results('Ipower');
%%
[x,y,z,ri,ci] = PM.identify_minimum('-Ipower');
P.LED_x = x;
P.Lens_ROC = y;
P.visible = 1;
P.beam_samples = 300;
if isfield(P,'rng_state')
  P = rmfield(P,'rng_state');
end
figure;
PM2 = WhiteCell(P);
%%
P.visible = 0;
P.beam_samples = 100;
P.evaluate_endpoints = 23; % 3 is M1, 23 is detector
P.rng_state = rng;
P.LED_x = -3;
PM = WhiteCell(P,'Lens_dx', .1:.1:1, 'Lens_ROC', 1.7:0.1:2.5);
figure;
PM.plot_results('Ipower');
%%
% Optimize Lens_dx for 24mm FL and LED_x
P.visible = 0;
P.beam_samples = 500;
P.evaluate_endpoints = 23; % 3 is M1, 23 is detector
P.rng_state = rng;
P.LED_x = -3;
PM = WhiteCell(P,'Lens_dx', .2:.05:.6, 'LED_x', -6:.1:-4.5);
figure;
PM.plot_results('Ipower');
