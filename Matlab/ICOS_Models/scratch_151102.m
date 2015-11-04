%%
% This model is to simulate work that Jordan is doing in the lab.
n = 2.4361;
P = ICOS_Model6.props;
P.R1 = 75;
P.R2 = 75;
P.r2 = P.r1;
P.y0 = 2.54;
P.dy = n*P.y0/P.R1;
P.dz = P.y0*sqrt(2)/P.R1;
P.Lenses{1} = 'Lens1';
P.Lens_Space = .2;
P.detector_spacing = 10;
P.Herriott_passes = 0;
P.herriott_spacing = 5;
P.HR = 0;
P.visibility = [0];
P.focus = 1;
P.visible = 1;
P.ICOS_passes_per_injection = 3;
%%
PM = ICOS_Model6(P);
%%
ddy = 0.1*linspace(-1,1,41);
ddz = linspace(-P.dz,.05,41);
P.visible = 0;
P.focus = 0;
P.ICOS_passes_per_injection = 60;
PM = ICOS_Model6(P,'dy',P.dy+ddy,'dz',P.dz+ddz);
%%
PM = PM.clean_results;
Y = (P.y0 - PM.Results.dy*P.mirror_spacing)/2.54;
Z = PM.Results.dz*P.mirror_spacing/2.54;
%%
mesh(Y,Z,double(PM.Results.eccentricity));
xlabel('Y in'); ylabel('Z in');
%%
contour(Y,Z,double(PM.Results.eccentricity),20);
xlabel('Y inches from axis'); ylabel('Z inches from axis');
set(gca,'XDir','Reverse');
grid;
title('Eccentricity vs Exit spot Location');
colorbar;
