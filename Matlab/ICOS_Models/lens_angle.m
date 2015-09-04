function [theta,ds] = lens_angle(r, d, s, f, x)
% r is the beam radius at x=0
% d is the beam divergence at x=0. Should be < 0
% s is the beam skew at x=0
% f is the focal length of the lens in question
% x is the position of the lens
% lens_angle is only relevant for converging focus, so d<0
% If the new divergence becomes positive, we've gone too far,
% so return a negative angle (past zero)
rx = sqrt((r+d.*x).^2 + s.^2.*x.^2);
dx = (d.*r + (d.^2+s.^2).*x)./rx - rx./f;
sx = s.*r./rx;
theta = rad2deg(atan(sqrt(sx.^2 + dx.^2)));
if dx > 0
  theta = -theta;
  if nargout > 1
    ds = 0;
  end
else
  if nargout > 1
    ds = -rx*dx/(dx^2+sx^2);
  end
end
