function draw_spots(R, r, phi, N)
% draw_spots(R, r, phi, N);
% Draw N spots at radius r phi radians apart.
% Also draw a circle of radius R for reference.
th = linspace(0, 2*pi, 100);
clf;
plot(R*sin(th),R*cos(th),'y');
hold on;
th = (0:N)*phi;
lap = th/(2*pi);
scatter(r*sin(th),r*cos(th),[],[1 1 1]);
hold off;
set(gca,'DataAspectRatio',[1 1 1],'color',[0,0,.3]);
% colorbar;
shg;
