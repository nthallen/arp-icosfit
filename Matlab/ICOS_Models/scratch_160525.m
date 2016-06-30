%%
r0 = 1.75*2.54/2;
RC = 0; % flat
P = 760; % Torr
E = 10.1e6; % psi
pr = 0.28;
CT = 0.5; % cm
O.C = [5 0 0];
O.D = [1 0 0];
R.O = [0 1 0];
dy = -0.01;
dz = 0.005;
R.D = [1 dy dz];
R.D = R.D/sqrt(sum(R.D.^2));
k = -(3/16)*P*(14.7/760)*(1-pr^2)/(E*CT^3*2.54);

t = 1;
X = @(t) R.O+t*R.D-O.C;
dxR = @(t) sum(X(t).*O.D);
rv = @(t) X(t)-dxR(t)*O.D;
r2 = @(t) sum(rv(t).^2);
if RC == 0
  dxr = @(t) k*(r0^2-r2(t))^2;
else
  dxr = @(t) sqrt(RC^2-r2(t))-RC + k*(r0^2-r2(t))^2;
end
ddxr = @(t) dxR(t)-dxr(t);
%%
t0 = fzero(ddxr,0);
