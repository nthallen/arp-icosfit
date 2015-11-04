%%
% Create models for Sayres' proposal
L = 35; % ICOS mirror spacing
B = 0.4; % Beam width, formerly known as 'W'
rd = .1; % detector radius (or half-width)
th = 15; % detector acceptance half-angle in degrees
C = 1000; % coherence length of laser
rdtanth = rd*tand(th);
mmin = ceil(C/(2*L));

k = 3;
R1 = 75;
%r1 = 1.596;
RL = 12.26;

k = 2;
R1 = 100;

m = mmin;
while(gcd(m,k) > 1)
  m = m + 1;
end
rmin = B/(2*sin(pi/m));
R2 = L*(R1-L)./(R1*sin(k*pi/m)^2-L);
r1_2 = rdtanth*sqrt(((R2-L).*R1.^2.*L)./((R1-L).*(R1+R2-L)));
r1max = sqrt(r1_2);
s1 = r1max.*sqrt((R1-L).*(R1+R2-L)./(R1.^2*L.*(R2-L)));
% SP.RL + SP.R1 + SP.R2 + SP.L + SP.Rw1
SP.R1 = R1;
SP.R2 = R2;
SP.L = L;
SP.Rw1 = B*.55;
SP.RL = SP.Rw1/s1;
Res = exparam(SP);
%%
Res.D1 = ceil((Res.r1+B)*4/2.54)/2; % round up to 1/2" diameter
Res.D2 = ceil((Res.r2+B)*4/2.54)/2; % round up to 1/2" diameter
Res.RD1 = ceil((Res.Rr1+B)*4/2.54)/2;
%%
P = render_model(Res);
P.n_overlap_spots = ceil(C/(2*L));
% P.Herriott_passes = 0;
P.ICOS_passes_per_injection = P.n_overlap_spots;
P.max_rays = 2000;

P.visible = 1;
P.focus = 1;
P.detector_spacing = 10;
PM = ICOS_Model6(P);
%%
IS = ICOS_search('mnc', 'sp1', 'RR1', Res.RR1, 'Rw1', Res.Rw1, ...
  'R1', Res.R1, 'L', Res.L, 'R2', Res.R2);
IS.search_ICOS_RIM;
%%
IS.search_focus2('select', 1, 'max_lenses', 1, 'focus_visible', 2);
%%
P = render_model(IS.res2(1),'max_rays',3000,'ICOS_passes_per_injection', ceil(C/(2*L)));
PM = ICOS_Model6(P);
%%
IB = ICOS_beam(@ICOS_Model6, P);
