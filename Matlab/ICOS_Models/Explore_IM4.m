%%
cd C:\Data\ES96\ICOS\3D
P = ICOS_Model4.props;

% Formatting
P.visible = false;
% P.avifile = 'L2_sweep.avi';
P.plot_endpoints = 0;
P.pause = 0;
P.title = 'Scan';
P.fmt_L2_X = 'L2\\_X = %.1f cm';
P.fmt_L1_R2 = 'L1\\_R2 = %.2f cm';

P.RC = 150;
P.herriott_spacing = 10; % Before the first ICOS mirror
P.HRC = 6*2.54; % Herriott radius of curvature
P.HR = 0;
P.injection_scale = 1;
P.y0 = P.Hr-0.54; % Location of Herriott hole
P.dy = .0344; % .04; replaced with optimal for circular alignment
P.dz = 0.0328; % 0.03; replace with optimal for circular alignment
P.dy = .0543; % optimal for RC = 150
P.dz = .0491; % optimal for RC = 150
% Optimized 4/22 further from optic edges
P.y0 = 1.905; % P.Hr/2
P.dy = 0.0308; % 
P.dz = 0.0286;
P.max_rays = 200;
P.ICOS_passes_per_injection = 70;
%%
P.visible = true;
P.visibility = [0 1 1 ];
P.injection_scale = 0.5;
PM = ICOS_Model4(P);
%PM = ICOS_Model4(P,'injection_scale', linspace(1,0.5,10));

%% This sequence is intended to look at what range of input angles
% will form successful ICOS alignments. Here I am assuming y0 is
% pretty close to the edge.
P.visible = false;
P.injection_scale = 1;
P.y0 = P.Hr/2;
PM = ICOS_Model4(P,'dy', linspace(-0.02,0.07,20),'dz',linspace(0,0.05,20));
PM = PM.clean_results;
PM.plot_results('eccentricity');
%%
P.visible = false;
[dy,dz,eccentricity] = PM.identify_minimum('eccentricity');
PM = ICOS_Model4(P,'dy', dy + linspace(-3.5e-3,+3.5e-3,20),'dz',dz + linspace(-2.5e-3,2.5e-3,20));
PM = PM.clean_results;
PM.plot_results('eccentricity');
[dy2,dz2,eccentricity2] = PM.identify_minimum('eccentricity');
%%
P.dy = dy2;
P.dz = dz2;
P.visible = true;
P.focus = 0;
PM = ICOS_Model4(P);

%%
% This is a parallel investigation, might want to skip
P.focus = 1;
P.dz = dz2;
PM = ICOS_Model4(P);
%%
[xyz,~,divergence,skew] = PM.M.extract_endpoints_skew(4);
x = 1:length(skew);
clf; plot(x,divergence,'.',x, skew,'*');
legend('divergence', 'skew');
%%
plot(xyz(:,2),xyz(:,3),'.r');
%%
r = sqrt(sum(xyz(:,[2 3]).^2,2));
plot(r,'.');
%%
fprintf(1,'mean(skew*r)/(mean(r)^2) = %.4f\n', ...
  mean(skew.*r)/mean(r.^2));
%%
% Now I'll reintroduce some eccentricity:
P.dz = P.dz/2;
PM = ICOS_Model4(P);
%%
[xyz,dxyz,divergence,skew] = PM.M.extract_endpoints_skew(4);
x = 1:length(skew);
r = sqrt(sum(xyz(:,[2 3]).^2,2));
clf; plot(r,divergence,'+',r, skew,'*');
legend('divergence', 'skew');
xlabel('r cm');
ylabel('cm');
fprintf(1,'mean(skew*r)/(mean(r)^2) = %.4f\n', ...
  mean(skew.*r)/mean(r.^2));

%% Now I'll try the same thing with a region farther from the rim
P.y0 = P.r/2;
PM = ICOS_Model4(P,'dy', linspace(-0.07,0.07,20),'dz',linspace(0,0.08,20));
PM = PM.clean_results;
PM.plot_results('eccentricity');
%%
[dy,dz,eccentricity] = PM.identify_minimum('eccentricity');
PM = ICOS_Model4(P,'dy', dy + linspace(-3.5e-3,+3.5e-3,20),'dz',dz + linspace(-2.5e-3,2.5e-3,20));
PM = PM.clean_results;
PM.plot_results('eccentricity');
[dy2,dz2,eccentricity2] = PM.identify_minimum('eccentricity');
%%
P.dy = dy2;
P.dz = dz2;
P.visible = true;
P.focus = 0;
PM = ICOS_Model4(P);

%% Now let's focus on the RIM. Let's start
cd C:\Data\ES96\ICOS\3D
P = ICOS_Model4.props;

% Formatting
P.visible = false;
% P.avifile = 'L2_sweep.avi';
P.plot_endpoints = 0;
P.pause = 0;
P.title = 'Scan';
P.fmt_L2_X = 'L2\\_X = %.1f cm';
P.fmt_L1_R2 = 'L1\\_R2 = %.2f cm';
P.fmt_herriott_spacing = 'RIM\\_X = %.2f cm';

P.herriott_spacing = 16.6; % Before the first ICOS mirror
P.RC = 150;
P.HRC = 24*2.54; % Herriott radius of curvature
P.HR = 1; % 0.98; % Herriott reflectivity
P.injection_scale = 1;
%P.y0 = P.Hr/2; % Location of Herriott hole
%P.dy = .0201; % circular injection for y0 = P.Hr/2
%P.dz = 0.0191; %
% Take these optimal values from section 1
%P.y0 = P.Hr-0.54; % Location of Herriott hole
%P.dy = .0543; % optimal for RC = 150
%P.dz = .0491; % optimal for RC = 150
% Optimized 4/22 further from optic edges
P.y0 = 1.905; % P.Hr/2
P.dy = 0.0308; % 
P.dz = 0.0286;
P.max_rays = 10000;
P.stop_ICOS = 1;
%%
P.visible = true;
P.stop_ICOS = 0;
P.max_rays = 600;
P.ICOS_passes_per_injection = 30;
PM = ICOS_Model4(P);
%%
P.visible = false;
P.injection_scale = 1;
P.plot_endpoints = 2;
PM = ICOS_Model4(P,'herriott_spacing',linspace(2,20,50));

%%
P.visible = false;
P.injection_scale = 0.5;
P.herriott_spacing = 10;
P.plot_endpoints = 3;
P.evaluate_endpoints = 3;
P.fmt_HRC = 'Herriott RoC = %.2f';
P.pause = 0;
P.stop_ICOS = 0;
P.ICOS_passes_per_injection = 63;
PM = ICOS_Model4(P,'HRC',[15.24 30.48 40.64 60.96 91.44 152.4],'herriott_spacing', 6:.2:20);
PM.clean_results;
%%
P.HRC = 91.44;
P.HRC = 15.24;
P.visible = true;
P.plot_endpoints = 0;
P.stop_ICOS = 0;
P.pause = 0;
P.ICOS_passes_per_injection = 63;
PM = ICOS_Model4(P,'HRC',[15.24 30.48 40.64 60.96 91.44 152.4]);
%%
PM.plot_results('mean_angle');
%%
PM.plot_results('max_radius');
%%
PM.plot_results('overlap');
%%
PM.plot_results('eccentricity');
%%
PM.plot_results('total_power');

%% Now let's try to focus
cd C:\Data\ES96\ICOS\3D
P = ICOS_Model4.props;

P.focus = 1;
P.ICOS_passes_per_injection = 63;
P.max_rays = 1000;
P.RC = 150;
P.L2_R1 = 2.123;
P.L2_R2 = 5.849;
P.L2_CT = 0.3;
P.L2_r = 1.27;
P.L2_X = 57.1;
% ZC-PM-25-25 1" positive meniscus d/f = 1
P.L3_r = 2.54/2;
P.L3_R1 = 2.291;
P.L3_R2 = 5.834;
P.L3_CT = 0.44;
P.L3_EFL = 2.54;
P.L3_X = P.L2_X + P.L2_CT + 0.4;
P.D_X = 60;

% Formatting
P.visible = true;
P.visibility = [0 0 0 1 1 1 1];
% P.avifile = 'L2_sweep.avi';
P.plot_endpoints = 0;
P.pause = 0;
P.title = 'Scan';
P.fmt_L2_X = 'L2\\_X = %.1f cm';
P.fmt_L1_R2 = 'L1\\_R2 = %.2f cm';
P.fmt_herriott_spacing = 'RIM\\_X = %.2f cm';

P.herriott_spacing = 16.6; % Before the first ICOS mirror
P.RC = 150;
P.HRC = 24*2.54; % Herriott radius of curvature
P.HR = 1; % 0.98; % Herriott reflectivity
P.injection_scale = 1;
%P.y0 = P.Hr/2; % Location of Herriott hole
%P.dy = .0201; % circular injection for y0 = P.Hr/2
%P.dz = 0.0191; %
% % Take these optimal values from section 1
% P.y0 = P.Hr-0.54; % Location of Herriott hole
% P.dy = .0543; % optimal for RC = 150
% P.dz = .0491; % optimal for RC = 150
% Optimized 4/22 further from optic edges
P.y0 = 1.905; % P.Hr/2
P.dy = 0.0308; % 
P.dz = 0.0286;
P.max_rays = 1000;

P.stop_ICOS = 0;
P.pause = 0;

% PM = ICOS_Model4(P, 'L2_X', linspace(60,55,15));
% PM = ICOS_Model4(P, 'D_X', linspace(59,61,20));
PM = ICOS_Model4(P);
%%
P.D_X = 61;
PM = ICOS_Model4(P);
%%
P.visible = false;
P.plot_endpoints = 7;
P.evaluate_endpoints = 7;
PM = ICOS_Model4(P, 'D_X', linspace(59,61,20));
%%
P.visible = false;
P.focus = 2;
P.visibility = [0 0 0 1 1 1 1 1 1];
P.L4_R1 = P.L2_R1;
P.L4_R2 = P.L2_R2;
P.L4_CT = P.L2_CT;
P.L4_X = 62;
P.L5_dX = 0.2; % arbitrary
P.D_X = P.L4_X + P.L4_CT + P.L5_dX + P.L5_CT + 2;
P.plot_endpoints = 8;
P.evaluate_endpoints = 8;
P.max_rays = 1000;
%%
PM = ICOS_Model4(P,'L4_X',linspace(62,59,10),'L5_dX',[0.2:.05:0.5]);
figure; PM.plot_results('max_radius');
%%
P.visible = true;
P.L4_X = 60;
P.L5_dX = 0.2;
P.D_X = P.L4_X + P.L4_CT + P.L5_dX + P.L5_CT + 2;
PM = ICOS_Model4(P);
%%
P.visible = false;
P.plot_endpoints = 9;
P.evaluate_endpoints = 9;
PM = ICOS_Model4(P,'L4_X', linspace(60,59,11));
%%
P.evaluate_endpoints = 9;
PM = ICOS_Model4(P,'L4_X',linspace(60,59,11),'L5_dX',[0.2:.05:0.5]);
figure; PM.plot_results('max_radius');
%%
P.visible = true;
P.L4_R1 = P.L2_R1;
P.L4_R2 = P.L2_R2;
P.L4_CT = P.L2_CT;
P.L4_X = 59;
P.L5_dX = 0.5;
PM = ICOS_Model4(P);

%%
P.visible = false;
P.evaluate_endpoints = 7;
P.max_rays = 1000;
PM = ICOS_Model4(P,'D_dY',[-0.3:.02:0.3]);
figure;
PM.plot_results('total_power'); shg;
%%
P.dy = 0.0201;
P.dz = 0.002;
P.mirror_spacing = 50;
P.max_rays = 1000;
P.ICOS_passes_per_injection = 10000;
PM = ICOS_Model4(P);
figure; PM.M.plot_endpoints(3);
%%
PM = ICOS_Model4(P,'D_dY',[-0.3:.02:0.3]);
figure; PM.plot_results('total_power'); shg;
%%
figure;
PM.M.plot_endpoints(7)


%% Look at overlap near 150/50
cd C:\Data\ES96\ICOS\3D
P = ICOS_Model4.props;

% Formatting
P.visible = false;
% P.avifile = 'L2_sweep.avi';
P.plot_endpoints = 0;
P.pause = 0;
P.title = 'Scan';
P.fmt_L2_X = 'L2\\_X = %.1f cm';
P.fmt_L1_R2 = 'L1\\_R2 = %.2f cm';

P.RC = 150;
P.herriott_spacing = 10; % Before the first ICOS mirror
P.HRC = 6*2.54; % Herriott radius of curvature
P.HR = 0;
P.injection_scale = 1;
P.y0 = P.Hr-0.54; % Location of Herriott hole
P.dy = .0344; % .04; replaced with optimal for circular alignment
P.dz = 0.0328; % 0.03; replace with optimal for circular alignment
P.dy = .0543; % optimal for RC = 150
P.dz = .0491; % optimal for RC = 150
P.max_rays = 200;
P.ICOS_passes_per_injection = 70;
P.focus = 0;
%%
P.evaluate_endpoints = 3;
P.visible = true;
PM = ICOS_Model4(P);
%%
P.injection_scale = 0.75;
P.visible = false;
PM = ICOS_Model4(P,'mirror_spacing', 48:.1:52);
PM.plot_results('overlap');


%% ------------------------------------------------------------------------------
% OK, new regime: lets go for eccentric
cd C:\Data\ES96\ICOS\3D
P = ICOS_Model4.props;
P.HR = 0;
P.visibility = [0 1 1];
P.visible = true;
P.injection_scale = 0.5;
P.focus = 1;
PM = ICOS_Model4(P);