%%
cd C:\Users\nort.ARP\Documents\SW\arp-icosfit\Matlab\ICOS_Models
%%
P = WhiteCell.props(25,2);
P.visible = 1;
P.beam_samples = 100;
P.beam_divergence = 2;
P.M0_Ddz = .2/25;
% P.M0_roc = 27;
P.Cell_Length = 25;
% P.visibility = [1 1 0 0];
P.evaluate_endpoints = 0;
PM = WhiteCell(P);
%xlim([-2 2]);
%%
% This analysis was moved into the model
rays = 1:PM.M.n_rays;
zs = find([ PM.M.Rays(rays).n_inc ] == 0);
pass = rays-interp1(zs,zs,rays,'previous','extrap')+1;
clear zs
%
% Now I want to calculate the spot size and power on each pass
max_pass = max(pass);
opt_n = zeros(max_pass,1);
opt_n(2:2:max_pass) = 2;
opt_n(4*P.N4) = 1;
opt_n(1:4:max_pass) = 3;
opt_n(3:4:max_pass) = 4;
%%
% This analysis was moved into evaluate_endpoints
n_opt = [PM.M.Rays(rays).n_opt];
Ipower = zeros(max_pass,1);
Opower = zeros(max_pass,1);
sizes = zeros(max_pass,1);
position = zeros(max_pass,3);
for i = 1:max_pass
  % fprintf(1,'Pass = %d\n', i);
  vf = find(pass == i);
  xyz = zeros(length(vf),3);
  for j=1:length(vf)
    vfi = vf(j);
    % fprintf(1, 'j = %d vf(j) = %d\n', j, vfi);
    ray = PM.M.Rays(vfi).ray;
    xyz(j,:) = ray.E;
    if opt_n(i) == n_opt(vfi) && ray.Inside
      Ipower(i) = Ipower(i) + ray.P;
    else
      Opower(i) = Opower(i) + ray.P;
    end
  end
  sizes(i) = sqrt(var(xyz(:,2))+var(xyz(:,3)));
  position(i,:) = mean(xyz) - PM.M.Optic{opt_n(i)}.O;
end
Ipower = Ipower/P.beam_samples;
Opower = Opower/P.beam_samples;
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
P.pause=0;
P.view = [-90,0];
PM = WhiteCell(P,'M2_fy', -(0:.01:0.5));
%%
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
PM = WhiteCell(P,'M12_Ddz',(-1:.04:1)*2);
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
%%
figure;
plot(power,'*');
%%
figure;
plot(spot_dia,'*r');
xlabel('Pass');
ylabel('Spot Diameter cm');
%%
figure;
plot(loss*100,'*r');
xlabel('Pass');
ylabel('% Loss');
