%%
% Testing equations for R1==R2 portion of the paper
L = 50;
r = 1;
W = .4;
rd = .1;
th = 15;
C = 3000;
%%
w1 = W/2;
m = 31;
k = 3;
rmin = W/(2*sin(pi/m));
phi = k*pi/m;
R = L/(1-cos(phi));
s = (rmin/R)*sqrt(2*R/L - 1);
%%
% SP.RL + SP.R1 + SP.R2 + SP.L + SP.Rw1
SP.R1 = R;
SP.R2 = R;
SP.L = L;
SP.Rw1 = W/2;
SP.RL = SP.Rw1/s;
%%
Res = exparam(SP);
%%
P = ICOS_Model6.props;
P.herriott_spacing = Res.RL;
P.HRC = Res.RR1;
P.Hr = Res.Rr1;
% P.HR = 0;
P.mirror_spacing = Res.L;
P.r1 = 2.54;
P.R1 = Res.R1;
P.R2 = Res.R2;
P.r2 = 2.54;
P.y0 = Res.r1;
P.dy = -Res.Rd2;
P.dz = Res.Rs2;
% P.Herriott_passes = 0;
P.ICOS_passes_per_injection = P.n_overlap_spots;
P.max_rays = 1000; % overkill, I know

P.focus = 0;
P.visibility = [0];
P.visible = 0;
%%
% PM = ICOS_Model6(P);
%%
IB = ICOS_beam(@ICOS_Model6, P);
%%
IB.Sample('opt_n', 3, 'n_optics', 3);
%%
Sample = IB.Res.NPass(1:IB.Res.N,:);
passes = 1:max(Sample(:,2));
spot = zeros(size(passes));
for i=1:length(passes)
  pass=passes(i);
  v = find(Sample(:,1)==1 & Sample(:,2) == pass);
  spot(i) = sqrt(sum(var(IB.Res.E(v,[2 3]))));
end
plot(passes,spot,'.'); shg;
