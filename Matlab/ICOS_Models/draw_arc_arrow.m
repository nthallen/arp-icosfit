function draw_arc_arrow(ax, O, r, phi, headlen, limrange)
% draw_arrow(ax, O, D, scale);
% ax Axes on which to draw
% O origin x,y
% r radius
% phi: 2-element vector indicating starting and ending angle
%  Angles are measured in radians clockwise from vertical
%  (by accident)
% headlen: arrowhead size.
% limrange: if specified, draws a limit line perpendicular to the arrowhead
if iscolumn(O)
  O = O';
end
if nargin < 6
  limrange = [];
end
th = linspace(phi(1), phi(2), 11);
X = O(1)+r*sin(th);
Y = O(2)+r*cos(th);
hold(ax,'on');
plot(ax,X,Y,'k');
draw_arrow(ax,[X(end-1) Y(end-1)], [X(end) Y(end)], headlen, limrange);
