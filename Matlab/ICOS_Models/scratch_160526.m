%% scratch_160526.m
% Explore mirror deflection using earlier solutions from HCl_TB4
load('IB_HCl_TB4_rd08_th13.0.1754_200_1000_102_23.5.16_50x100.mat');
load('SR_HCl_TB4_rd08_th13.0.mat');
SRindex = 1754; % could be extracted from IB mnc
dBL = SR.Summary(1754).dBL;
%%
% Start with basic tests of ICOS
P = IB.P;
P.focus = 0;
P.HR = 0;
P.visible = 0;
P.CT1 = 0.5;
P.CT2 = 0.5;

%%
delta = 0.002;
%%
while delta > 1e-5
  ddy = delta*linspace(-1,1,21);
  ddz = delta*linspace(-1,1,21);
  P.visible = 0;
  P.focus = 0;
  P.evaluate_endpoints = 3;
  P.ICOS_passes_per_injection = 60;
  PM = ICOS_Model6(P,'dy',P.dy+ddy,'dz',P.dz+ddz);
  %
  PM = PM.clean_results();
  PM.plot_results('eccentricity');
  figure;
  PM.plot_results('overlap');
  %
  [ new_dy, new_dz ] = PM.identify_minimum('eccentricity');
  % fprintf(1,' ddY = %f   ddZ = %f\n', new_dy - P.dy, new_dz - P.dz);
  fprintf(1,'TddY = %f  TddZ = %f\n', new_dy - IB.P.dy, new_dz - IB.P.dz);
  %%
  P.dy = new_dy;
  P.dz = new_dz;
  delta = delta/5;
end
% This is now the baseline configuration:
%   Original cavity length, corrected mirror thickness, optimized dy,dz
Pbl = P;
%%
Pbl.evaluate_endpoints = 2;
Pbl.plot_endpoints = 0;
dLv = linspace(2*dBL(1),2*dBL(2),61);
PM = ICOS_Model6(Pbl,'mirror_spacing', Pbl.mirror_spacing + dLv);
% is the overlap==0 region contiguous? Does it agree with dBL?
iszero = PM.Results.overlap == 0;
if iszero(1) == 0 && iszero(end) == 0 && sum(diff(iszero) ~= 0) == 2
  fprintf(1,'overlap==0 region is contiguous\n');
else
  fprintf(1,'overlap==0 region is NOT as expected\n');
end
%%
edBL = minmax(PM.Results.mirror_spacing(iszero))-Pbl.mirror_spacing;
ddBL = diff(dBL);
eddBL = diff(edBL);
fprintf(1,'eddBL differs by %.1f%%\n', 100*(eddBL-ddBL)/ddBL);
for i=1:2
  fprintf(1,'edBL(%d) differs by %.1f%%\n', i, 100*(edBL(i)-dBL(i))/dBL(i));
end
%
% P.plot_endpoints = 3;
% P.ICOS_passes_per_injection = 200;
% PM = ICOS_Model6(P);
%
% figure;
% P.dP = 0;
% P.visible = 0;
% PM = ICOS_Model6(P);
%
% figure;
% P.dP = 760;
% PM = ICOS_Model6(P);
%%
P.dP = 760;
P.evaluate_endpoints = 2;
dLv = linspace(2*dBL(1),2*dBL(2),61);
PM = ICOS_Model6(P,'mirror_spacing', P.mirror_spacing + dLv);
% is the overlap==0 region contiguous? Does it agree with dBL?
iszeroP = PM.Results.overlap == 0;
if iszeroP(1) == 0 && iszeroP(end) == 0 && sum(diff(iszeroP) ~= 0) == 2
  fprintf(1,'overlapP==0 region is contiguous\n');
else
  fprintf(1,'overlapP==0 region is NOT as expected\n');
end
%%
figure;
plot(PM.Results.mirror_spacing(iszero)-P.mirror_spacing,PM.Results.overlap(iszero),'.r', ...
  PM.Results.mirror_spacing(~iszero)-P.mirror_spacing,PM.Results.overlap(~iszero),'.b');
%%
figure;
plot(PM.Results.mirror_spacing(iszero)-P.mirror_spacing,PM.Results.eccentricity(iszero),'.r', ...
  PM.Results.mirror_spacing(~iszero)-P.mirror_spacing,PM.Results.eccentricity(~iszero),'.b');
%%
% Now select the best eccentricity that is more than 10/1000" from either
% overlap region
zor = minmax(PM.Results.mirror_spacing(iszero&iszeroP))-P.mirror_spacing;
zor = zor + .01*2.54*[1,-1];
dL = PM.Results.mirror_spacing-P.mirror_spacing;
zorOK = dL >= zor(1) & dL <= zor(2);
zorE = PM.Results.eccentricity(zorOK);
minEzn = find(zorE == min(zorE));
minEn = find(zorOK);
minEn = minEn(minEzn);
optL = PM.Results.mirror_spacing(minEn);
%%
figure;
plot(PM.Results.mirror_spacing(iszero)-P.mirror_spacing, ...
  PM.Results.eccentricity(iszero),'.r', ...
  PM.Results.mirror_spacing(~iszero)-P.mirror_spacing, ...
  PM.Results.eccentricity(~iszero),'.b', ...
  optL-P.mirror_spacing,zorE(minEzn),'*g');
%%
% Now we've selected a spot in the overlap region
% Update mirror_spacing, then reoptimize dy,dz for eccentricity
P.mirror_spacing = optL;
delta = 0.002;
%
while delta > 1e-5
  ddy = delta*linspace(-1,1,21);
  ddz = delta*linspace(-1,1,21);
  P.visible = 0;
  P.focus = 0;
  P.evaluate_endpoints = 3;
  P.ICOS_passes_per_injection = 60;
  PM = ICOS_Model6(P,'dy',P.dy+ddy,'dz',P.dz+ddz);
  %
  PM = PM.clean_results();
  PM.plot_results('eccentricity');
  figure;
  PM.plot_results('overlap');
  %
  [ new_dy, new_dz ] = PM.identify_minimum('eccentricity');
  fprintf(1,' ddY = %f   ddZ = %f\n', new_dy - P.dy, new_dz - P.dz);
  fprintf(1,'TddY = %f  TddZ = %f\n', new_dy - Pbl.dy, new_dz - Pbl.dz);
  %%
  P.dy = new_dy;
  P.dz = new_dz;
  delta = delta/5;
end
% Now we have the ICOS-optimized configuration for 760 Torr
fprintf(1,'Baseline 760 Difference\n');
fprintf(1,'L: %7.2f %7.2f %7.2f\n', Pbl.mirror_spacing, P.mirror_spacing, ...
  P.mirror_spacing - Pbl.mirror_spacing);
fprintf(1,'dy: %7.5f %7.5f %7.5f\n', Pbl.dy, P.dy, ...
  P.dy - Pbl.dy);
fprintf(1,'dz: %7.5f %7.5f %7.5f\n', Pbl.dz, P.dz, ...
  P.dz - Pbl.dz);
%%
save s160526.5.mat Pbl P dBL
%%
load s160526.5.mat
%%
% Now run ICOS_beam analysis on Pbl, P w/ dP=[] and P w/ dP=760,
% all using the same seed.
Pbl.evaluate_endpoints = 0;
Pbl.focus = 1;
IBopt = {'beam_samples',100,'ICOS_passes',50,'n_optics',6, ...
  'opt_n', 6, 'Track_Power', 1};
IBbl = ICOS_beam(@ICOS_Model6, Pbl);
IBbl.Sample('mnc', 'dP_bl.5', IBopt{:});
IBbl.Integrate;
IBbl.savefile;
%%
P.dP = [];
P.evaluate_endpoints = 0;
P.focus = 1;
IBP0 = ICOS_beam(@ICOS_Model6, P);
IBP0.Sample('mnc', 'dP_P0.5', 'rng_state', IBbl.IBP.rng_state, IBopt{:});
IBP0.Integrate;
IBP0.savefile;
%%
P.dP = 760;
IBP760 = ICOS_beam(@ICOS_Model6, P);
IBP760.Sample('mnc', 'dP_P760.5', 'rng_state', IBbl.IBP.rng_state, IBopt{:});
IBP760.Integrate;
IBP760.savefile;
