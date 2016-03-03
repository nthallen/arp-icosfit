function animate_interleave(m,k)
% figure_overlap
% Draws a figure to explain the overlap calculation
n = m;
B = .4;
w = B/2;
phi = pi/n;
Phi = k*pi/n;
rmin = w/sin(phi);
rmax = 2*rmin;
N = 1:n;
% phis = 2*N*phi; % angle of spot centers, starting at vertical
% O = rmin*[sin(phis') cos(phis')];
% for i=1:n
%   draw_circle(O(i,:), w, 'g');
% end

dN = evaluate_interleave(m,k);
phi_min = asin(B./(2*rmax));
Phi_n = Phi + (phi_min-phi)/dN(2);
Phi_p = Phi + (phi_min-phi)/dN(1);

f1 = figure;
phis3 = 2*N*Phi; % same for outer spots
O3 = rmax*[sin(phis3') cos(phis3')];
for i=1:n
  draw_circle(O3(i,:), w, 'b', i);
end
set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
hold off;

f2 = figure;
phis3 = 2*N*Phi_n; % same for outer spots
O3 = rmax*[sin(phis3') cos(phis3')];
for i=1:n
  draw_circle(O3(i,:), w, 'b', i);
end
set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
hold off;

f3 = figure;
phis3 = 2*N*Phi_p; % same for outer spots
O3 = rmax*[sin(phis3') cos(phis3')];
for i=1:n
  draw_circle(O3(i,:), w, 'b', i);
end
set(gca,'DataAspectRatio',[1 1 1],'XTick',[],'YTick',[]);
hold off;


function draw_circle(O, r, c, n)
  sth = linspace(0,2*pi,21);
  x = O(1)+r*sin(sth);
  y = O(2)+r*cos(sth);
  plot(x,y,c);
  h = text(O(1),O(2),sprintf('%d',n));
  set(h,'HorizontalAlignment','center','VerticalAlignment','middle');
  hold on;
