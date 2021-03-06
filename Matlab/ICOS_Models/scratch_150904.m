%%
% Interleave 2, r1 == r2
% Look at engineering an ideal cell without restriction to existing parts
% This version is proposing an interleaved spot pattern in the ICOS cell
% sropt_e: Trying r1 > sqrt(srL/sin(phi))
mnc = 'sropt_e';
W = .4; % 4 mm beam diameter
C = 3000; % 30 m coherence length
L = 50; % 50 cm cell length
r_d = 0.08; % detector radius 1 mm, reduced for optimal fill
th_d = 14.9; % detector acceptance angle 15 degrees
sr = r_d*tand(th_d);
m = ceil(C/(4*L)-1/2);
phi = 2*pi/(2*m+1);
r1 = sqrt(sr*L/sin(phi))*1.2; % actually r1 > sqrt(...)
w1 = r1*sin(phi);
s1 = sr/r1;
w2 = s1*L;
D1 = (r1+W)*2/2.54; % Diameter in inches
h1 = sqrt(r1^2-w1^2);
%%
Rs2 = s1;
Rw1 = W*.5;
RL = Rw1/Rs2;
%%
w2 = s1*L;
r2 = r1*w2/w1;
h2 = h1*r2/r1;
R1 = L*r1/(r1-h2);
R2 = L*r2/(r2-h1);
%%
% Now I've got R1, R2, L, r1, RL (plus sr)
% Also w1 and Rw1, s1, Rs2
n = 2.4361;
d1 = r1/R1;
RR2 = -R1/n;
Rd2 = -n*d1;
Rr2 = r1;
Rh1 = Rr2*(RR2-RL)/RR2;
Rr1 = sqrt(Rw1^2+Rh1^2);
Rh2 = Rh1*Rr2/Rr1;
RR1 = Rr1*RL/(Rr1-Rh2);
%%
% or R1, R2, RR1, L, Rw1
IS = ICOS_search('mnc', mnc, 'R1', R1, 'R2', R2, 'RR1', RR1, ...
  'L', L, 'Rw1', Rw1, 'RD1_margin', 5);
IS.search_ICOS_RIM;
%%
IS.search_focus2('select',1,'focus_visible',2,'max_lenses',2);
openvar('IS');
%%
IS.analyze;
%%
P = render_model(IS.res2(1));
PM = ICOS_Model6(P);
% clear P;
% P.R1 = R1;
% P.R2 = R2;
% P.RR1 = RR1;
% P.L = L;
% P.Rw1 = Rw1;
% res = exparam(P);
%%
mnc = strrep(IS.ISopt.mnc,'_','\_');

for i=1:5
  render_model(IS.res2(i),'view',[0 0],'ICOS_passes_per_injection',31,'max_rays',20000);
  title(sprintf('%s: Nres2: %d', mnc, IS.res2(i).Nres2));
end
%%
NM = zeros(length(IS.res2),1);
for i=1:length(NM)
  re = regexp(IS.res2(i).Lenses, '_NM_');
  re = [ re{:} ];
  if ~isempty(re)
    NM(i) = 1;
  end
end

