function draw_spots(R, r, phi, N, B)
% draw_spots(R, r, phi, N, B);
% Draw N spots at radius r phi radians apart.
% Also draw a circle of radius R for reference.
% If present use B as the spot diameter
th = linspace(0, 2*pi, 100);
clf;
plot(R*sin(th),R*cos(th),'y');
hold on;
th = (0:N)*phi;
lap = th/(2*pi);
if nargin < 5
  scatter(r*sin(th),r*cos(th),[],[1 1 1]);
else
  spth = linspace(0, 2*pi, 21);
  spx = B*sin(spth)/2;
  spy = B*cos(spth)/2;
  for i=1:length(th)
    plot(r*sin(th(i))+spx, r*cos(th(i))+spy, 'w');
  end
end
hold off;
set(gca,'DataAspectRatio',[1 1 1],'color',[0,0,.3]);
% colorbar;
shg;
