function [x,theta_o,ds_o] = pick_lens_x2(r, d, s, f, theta, lens_r)
% r is the beam radius at x=0
% d is the beam divergence at x=0. Should be < 0
% s is the beam skew at x=0
% f is the focal length of the lens in question
% theta is the target max angle
% lens_r is the lens radius, but the caller should subtract any desired
%   margin.
% Returns:
% x: empty if no solution. Scalar if only a boundary solution
%    [x0 x1] if only a boundary condition. In this case x0
%      is the closest to the desired focus, but values between
%      x0 and x1 are possible.
%    [x0 xmin xmax] if a good solution is found. x0 is pretty close
%      to the desired focus and lies between xmin and xmax.
%-----
% Add a lens radius parameter, lens_r. If r > lens_r, set xmin
% to the point where rx < lens_r-.3
x = [];
if nargout > 1
  theta_o = [];
  ds_o = [];
end
rmax = lens_r;
if d >= 0 && (f < 0 || r >= rmax)
  return; % no solutions with this lens. Heading the wrong direction
    % or Too small
end
V = [
  d^2+s^2
  2*r*d
  r^2-rmax^2
  ];
rts = roots(V);
rts = rts(imag(rts)==0);
rts = rts(rts > 0);
rts = min(rts); % rts should be where the radius == rmax
if isempty(rts)
  fprintf(1,'pick_lens_x2(%f,%f,%f,%f) no roots\n', r, d, s, rmax);
  return;
end
if d >= 0
  xmin = 0.1; % Since we can assume f>0, xmin represents the minimal correction
  xmax = rts; % and xmax will be the maximal correction.
else
  th0 = atand(sqrt(d^2+s^2));
  if ( th0 < theta && f < 0 ) || (th0 > theta && f > 0)
    return; % Heading the wrong direction
  end
  if r > rmax
    rmin = r*s/sqrt(s^2+d^2);
    if rmin < rmax
      xmin = rts;
    else
      return;
    end
  else
    xmin = 0.1;
  end
  if f > 0
    xmax = -r*d/(d^2 + s^2);
  else
    V = [
      -(d^2+s^2)^2/f
      -3*r*d*(d^2+s^2)/f
      -r^2*(3*d^2+s^2)/f
      r^2*(s^2-r*d/f)
      ];
    nfrts = roots(V);
    nfrts = nfrts(imag(nfrts)==0);
    nfrts = nfrts(nfrts > 0);
    nfrts = min(nfrts); % rts is where dx-rx/f is minimum
    if isempty(nfrts)
      error('MATLAB:HUARP:noroot','Expected root for f<0');
    end
    xmax = nfrts;
  end
end

try
  x0 = fzero(@(X) lens_angle(r, d, s, f, X)-theta, [xmin xmax]);
  x = [x0 xmin xmax];
  if nargout > 1
    [theta_o,ds_o] = lens_angle(r,d,s,f,x0);
  end
catch
  % Did not reach theta. Cases:
  % d > 0 => f > 0
  %   xmin: smallest correction
  %   xmax: largest correction
  %   If both overshoot (la0 > theta), take xmin
  %   else take xmax (both undershoot)
  % d < 0 && th0 < theta => f > 0
  %   xmin: largest correction
  %   xmax: smallest correction
  %   If both overshoot (la0 > theta), take xmax
  %   else take xmin
  % d < 0 && th0 > theta => f < 0
  %   xmin: largest correction (could go to d>0, la0 < 0)
  %   xmax: smallest correction
  %   if both overshoot (la0 < theta), take xmax
  %   else take xmin
  x = [xmin xmax];
  [la0,ds] = lens_angle(r,d,s,f,xmin);
  if d > 0
    if la0 < theta
      [la0,ds] = lens_angle(r,d,s,f,xmax);
      x = [xmax xmin];
    end
  elseif f > 0
    if la0 > theta
      [la0,ds] = lens_angle(r,d,s,f,xmax);
      x = [xmax xmin];
    end
  elseif la0 < theta
    [la0,ds] = lens_angle(r,d,s,f,xmax);
    x = [xmax xmin];
  end
  if nargout > 1
    theta_o = la0;
    ds_o = ds;
  end
end
