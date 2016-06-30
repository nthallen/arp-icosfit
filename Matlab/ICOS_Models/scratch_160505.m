%% Mirror displacement
R = 200/2.54; % Radius of curvature in inches
r0 = 1.75/2; % Mirror radius in inches
P = (760-150)*14.7/760; % Pressure differential in psi
E = 10.1e6; % Young's modules, psi
t = 0.5/2.54; % Mirror thickness (6mm) in inches
v = 0.28; % Poisson's ratio
m = 1/v;
%%
r = (0:.01:1)'*r0;
k = -3*P*(m^2-1)/(16*E*m^2*t^3);
dx = k*(r0^2-r.^2).^2;
x0 = sqrt(R^2-r.^2) - R; 
x = x0 + dx;
%%
plot(x0,r,x,r);
xlabel('in');
ylabel('in');
%%
plot(dx*2.54,r);
xlabel('cm');
ylabel('in');
%%
% The slope of the ray normal to the mirror surface
% is the negative of the slope of mirror surface when we swap
% coordinates. (It would be the negative reciprocal otherwise)
dxdr = -r./sqrt(R^2-r.^2)-4*r.*k.*(r0^2-r.^2);
drdx = -dxdr;
R1 = -x + r./drdx;
%%
plot(r,R1,r,ones(size(R1))*R);
%%
plot(r,100*(R1-R)/R);
ylabel('Percent Change');
xlabel('radius in');
