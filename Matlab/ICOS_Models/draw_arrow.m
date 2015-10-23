function draw_arrow(ax, O, D, headlen, limrange)
% draw_arrow(ax, O, D, scale);
% ax Axes on which to draw
% O origin x,y
% D destination x,y (where the arrowhead goes)
% headlen: arrowhead size. Defaults to 0.1 * the length
% limrange: if specified, draws a limit line perpendicular to the arrowhead
if iscolumn(O)
  O = O';
end
if iscolumn(D)
  D = D';
end
len = sqrt(sum((D-O).^2));
if nargin < 4
  headlen = len*0.1;
end
X = (O-D)/len;
Y = X * [0 1; -1 0];
if Y(2) < 0
  Y = -Y; % Make sure the 'Y' direction is upwards
end
x = [len 0 headlen headlen*0.8 headlen 0];
y = [0   0 headlen*0.2 0 -headlen*0.2 0];
P = [x' y'];
M = [X; Y];
Pr = P*M + ones(length(x),1)*D;
hold(ax,'on');
plot(ax, Pr(:,1), Pr(:,2), 'k');
if nargin >= 5 && ~isempty(limrange)
  if iscolumn(limrange)
    limrange = limrange';
  end
  x = [0 0];
  y = limrange;
  P = [x' y'];
  Pr = P*M + ones(length(x),1)*D;
  plot(ax, Pr(:,1), Pr(:,2), 'k');
end
hold(ax,'off');
