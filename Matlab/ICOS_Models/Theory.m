%% Theory
% R = radius of curvature for ICOS mirrors
% r = axial radius of spot pattern on mirror for circular alignment
% d = mirror spacing at axis
% theta = angle for one pass
R = 150;
r = 1.5 * 2.54 - 0.3;
d = 50;

% For a circular alignment, P0 and P2 are on the left mirror (M1),
% P1 is on the right mirror (M2)
% M1 is centered at (0,0,0), center of curvature at (R,0,0)
% M2 is centered at (d,0,0), center of curvature at (d-R,0,0)
% Circular pattern is offset in x by R-sqrt(R^2-r^2) from the mirror origin
% P1 is at (d+sqrt(R^2-r^2)-R, 0, r)
% P0 is at (R-sqrt(R^2-r^2), sqrt(r^2-h^2), h)
% h = z coordinate of P0
h = r*(1 - d/sqrt(R^2 - r^2));
P0 = [R-sqrt(R^2-r^2), sqrt(r^2-h^2), h];
P1 = [d+sqrt(R^2-r^2)-R, 0, r];
%%
D = P1-P0;
D = D/D(1);

div = sum(D([2 3]).*P1([2 3]))/r;
skew = sum(D([2 3]).*[-P1(3) P1(2)])/r;

%%
% Now let's do simple focus math. If we have a specific divergence and
% skew at P1, we have:
%   [y,z,dy,dz] = [0,r,skew,div]
% After a focusing element of focal length f, we get:
%   [y,z,dy,dz] = [0,r,skew,div-r/f]
% After propagating a distance x, that becomes:
%   [y,z,dy,dz] = [sx,r+x(d-r/f),s,d-r/f]
% Minimum radius is achieved when x = [r*(r/f-d)]/[s^2 + (r/f-d)^2]
% Minimum radius is:
%   minrad = s*r/sqrt(s^2+(r/f-d)^2)
dl = 0.1; % radius of detector
fmax = r/(div + (skew/dl)*sqrt(r^2-dl^2));
