%%
% Look at engineering an ideal cell without restriction to existing parts
W = .4; % 4 mm beam diameter
C = 3000; % 30 m coherence length
L = 50; % 50 cm cell length
r_d = 0.09; % detector radius 1 mm
th_d = 14.9; % detector acceptance angle 15 degrees
w1 = W/2;
r1 = W*C/(4*pi*L);
D1 = (r1+W)*2/2.54; % Diameter in inches
sr = r_d*tand(th_d);
s1 = sr/r1;
h1 = sqrt(r1^2-w1^2);
%%
Rs2 = s1;
Rw1 = W*.7;
RL = Rw1/Rs2;
%%
% Suppose R1=75. R2 = -140.9, pretty close to ZC-PX-25-100 at -139.96
R1 = 75;
k = r1^4*(R1-L)/(R1^2*L*sr^2);
R2 = (k/(1-k))*R1+L;
%%
% Now I've got R1, R2, L, r1, RL (plus sr)
% Also w1 and Rw1, s1, Rs2
n = 2.4361;
r2 = r1*sqrt(R2*(R1-L)/(R1*(R2-L)));
%h1 = r1*(R1-L)*(R2-L)/(R1*R2);
h2 = h1*r2/r1;
% w1 = sqrt(r1^2-h1^2);
w2 = w1*r2/r1;
% s1 = w2/L;
d1 = r1/R1;
RR2 = -R1/n;
Rd2 = -n*d1;
Rr2 = r1;
% Rs2 = s1;
% Rw1 = Rs2*RL;
Rh1 = Rr2*(RR2-RL)/RR2;
Rr1 = sqrt(Rw1^2+Rh1^2);
Rh2 = Rh1*Rr2/Rr1;
RR1 = Rr1*RL/(Rr1-Rh2);
RD1 = Rr1*2/2.54;
%%
% or R1, R2, RR1, L, Rw1
IS = ICOS_search('mnc', 'sropt4', 'R1', R1, 'R2', R2, 'RR1', RR1, 'L', L, 'Rw1', Rw1, 'RD1_margin', 5);
IS.search_ICOS_RIM;
%%
IS.search_focus;
%%
IS.analyze('select',2);
%%
P = render_model(IS.res2(2));
PM = ICOS_Model6(P);
% clear P;
% P.R1 = R1;
% P.R2 = R2;
% P.RR1 = RR1;
% P.L = L;
% P.Rw1 = Rw1;
% res = exparam(P);


