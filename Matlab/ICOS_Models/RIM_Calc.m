%% Check out optimal RIM
R0 = 75;
r0 = 3.5;
d0 = r0/R0;
L0 = 50;
s0 = (r0/R0)*sqrt(2*R0/L0 - 1);
n = 2.4361;
R1 = -R0/n;
r1 = r0;
d1 = d0*n;
s1 = s0;
L1 = .1:.1:60;
e2 = L1*s1;
h2 = r1 + d1*L1;
r2 = sqrt(e2.^2 + h2.^2);
h1 = r1.*h2./r2;
R2 = (r2.*L1)./(r2-h1);
figure; plot(L1,R2,'*');

R = [30.48 40.64 60.96];
L = interp1(R2,L1,R);
hold on;
plot(L,R,'*r');
hold off;

%% This calculation is for when the RIM's radius of curvature is known
% (i.e. selected from a catalog)
% r2, R2, d2 and s2 are determined from the ICOS configuration
% These parameters are for IM5 75x75 cylindrical orientation
r2 = 2;
R2 = -30.8;
d2 = -0.065;
s2 = 0.038;
FL1 = 8;

R1 = FL1*2*2.54;
sd2 = s2^2+d2^2;
a = sd2;
b = -2*r2*d2-R1*sd2;
c = r2^2+R1*r2*d2;
disc = b^2-4*a*c;
L = (-b+sqrt(disc))/(2*a);
L2 = (-b-sqrt(disc))/(2*a);
w1 = s2*L;
h1 = r2-d2*L;
th = rad2deg(atan(2*w1/h1));
Hnp = 360/th;
r1 = sqrt(h1^2+w1^2);

%% IM6r2 RIM
% r2, R2, d2 and s2 are determined from the ICOS configuration
% These parameters are for IM6r2 75x-28.12 conical orientation
r2 = 2;
R2 = -30.8;
d2 = -0.065;
s2 = 0.00377;
% Ideal FL is 16.5 for 2 = 0.2, diameter above 4"
% Edmund lists: Dia x FL in inches as:
%  5x4, 6x12, 6x24, 8x24
FL1 = 24;

R1 = FL1*2*2.54;
sd2 = s2^2+d2^2;
a = sd2;
b = -2*r2*d2-R1*sd2;
c = r2^2+R1*r2*d2;
disc = b^2-4*a*c;
L = (-b+sqrt(disc))/(2*a);
L2 = (-b-sqrt(disc))/(2*a);
w1 = s2*L;
h1 = r2-d2*L;
th = rad2deg(atan(2*w1/h1));
Hnp = 360/th;
r1 = sqrt(h1^2+w1^2);

%% IM6r3 RIM
% r2, R2, d2 and s2 are determined from the ICOS configuration
% These parameters are for IM6r2 75x-28.12 conical orientation
r2 = 3;
R2 = -30.8;
d2 = -0.097;
s2 = 0.00565;
% Ideal FL is 16.5 for 2 = 0.2, diameter above 4"
% Edmund lists: Dia x FL in inches as:
%  5x4, 6x12, 6x24, 8x24
FL1 = 12;

R1 = FL1*2*2.54;
sd2 = s2^2+d2^2;
a = sd2;
b = -2*r2*d2-R1*sd2;
c = r2^2+R1*r2*d2;
disc = b^2-4*a*c;
L = (-b+sqrt(disc))/(2*a);
L2 = (-b-sqrt(disc))/(2*a);
w1 = s2*L;
h1 = r2-d2*L;
th = rad2deg(atan(2*w1/h1));
Hnp = 360/th;
r1 = sqrt(h1^2+w1^2);

%% IM6ar2 RIM
% r2, R2, d2 and s2 are determined from the ICOS configuration
% These parameters are for IM6r2 75x-28.12 conical orientation
r2 = 2;
R2 = -30.8;
d2 = -0.065;
s2 = 0.0115;
% Ideal FL is 16.5 for 2 = 0.2, diameter above 4"
% Edmund lists: Dia x FL in inches as:
%  3x8, 3x12, 5x4, 6x12, 6x24, 8x24
FL1 = 10;

R1 = FL1*2*2.54;
sd2 = s2^2+d2^2;
a = sd2;
b = -2*r2*d2-R1*sd2;
c = r2^2+R1*r2*d2;
disc = b^2-4*a*c;
L = (-b+sqrt(disc))/(2*a);
L2 = (-b-sqrt(disc))/(2*a);
w1 = s2*L;
h1 = r2-d2*L;
th = rad2deg(atan(2*w1/h1));
Hnp = 360/th;
r1 = sqrt(h1^2+w1^2);

%% IM6br2 RIM
% r2, R2, d2 and s2 are determined from the ICOS configuration
% These parameters are for IM6r2 75x-28.12 conical orientation
r2 = 2;
R2 = -30.8;
d2 = -0.065;
s2 = 0.0147;
% Ideal FL is 16.5 for 2 = 0.2, diameter above 4"
% Edmund lists: Dia x FL in inches as:
%  3x8, 3x12, 5x4, 6x12, 6x24, 8x24
FL1 = 12;

R1 = FL1*2*2.54;
sd2 = s2^2+d2^2;
a = sd2;
b = -2*r2*d2-R1*sd2;
c = r2^2+R1*r2*d2;
disc = b^2-4*a*c;
L = (-b+sqrt(disc))/(2*a);
L2 = (-b-sqrt(disc))/(2*a);
w1 = s2*L;
h1 = r2-d2*L;
th = rad2deg(atan(2*w1/h1));
Hnp = 360/th;
r1 = sqrt(h1^2+w1^2);
