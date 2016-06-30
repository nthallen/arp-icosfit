function H = draw_arrow3(ax, O, D, headlen, limrange)
% H = draw_arrow(ax, O, D, headlen[, limrange]);
% ax Axes on which to draw
% O origin x,y
% D destination x,y (where the arrowhead goes)
% headlen: arrowhead size. Defaults to 0.1 * the length
% limrange: if specified, draws a limit line perpendicular to the arrowhead
% Returns handle to the line objects created
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
[az,el] = view;
vw = [cosd(el)*sind(az),-cosd(el)*cosd(az),sind(el)];
Y = cross(X,vw);
x = [len 0 headlen headlen*0.8 headlen 0];
y = [0   0 headlen*0.2 0 -headlen*0.2 0];
z = [0 0 0];
P = [x' y'];
M = [X; Y];
Pr = P*M + ones(length(x),1)*D;
hold(ax,'on');
h1 = plot3(ax, Pr(:,1), Pr(:,2), Pr(:,3), 'k');
if nargin >= 5 && ~isempty(limrange)
  if iscolumn(limrange)
    limrange = limrange';
  end
  x = [0 0];
  y = limrange;
  P = [x' y'];
  Pr = P*M + ones(length(x),1)*D;
  h2 = plot(ax, Pr(:,1), Pr(:,2), 'k');
  if nargout > 0
    H = [h1; h2];
  end
else
  if nargout > 0
    H = h1;
  end
end
hold(ax,'off');
