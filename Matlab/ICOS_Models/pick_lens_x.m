function [x,theta_o,ds_o] = pick_lens_x(r, d, s, f, theta, lens_r)
% r is the beam radius at x=0
% d is the beam divergence at x=0. Should be < 0
% s is the beam skew at x=0
% f is the focal length of the lens in question
% theta is the target max angle
%-----
% Add a lens radius parameter, lens_r. If r > lens_r, set xmin
% to the point where rx < lens_r-.3
x = [];
if nargout > 1
  theta_o = [];
  ds_o = [];
end
if d < 0
  lens_r = lens_r - 0.3;
  xmax = -r*d/(d^2 + s^2);
  if r > lens_r
    rmin = r*s/sqrt(s^2+d^2);
    if rmin < lens_r
      V = [
        d^2+s^2
        2*r*d
        r^2-lens_r^2
        ];
      rts = roots(V);
      rts = rts(imag(rts)==0);
      rts = rts(rts > 0);
      rts = min(rts);
      if isempty(rts)
        error('MATLAB:HUARP:ExpectedRoot','Expected to find a root');
      end
      xmin = rts;
    else
      return;
    end
  else
    xmin = 0.1;
  end
  try
    x = fzero(@(X) lens_angle(r, d, s, f, X)-theta, [xmin xmax]);
    if nargout > 1
      [theta_o,ds_o] = lens_angle(r,d,s,f,x);
    end
  catch
    [la0,ds] = lens_angle(r,d,s,f,xmin);
    if la0 >= 0 && la0 <= theta
      x = xmin;
      if nargout > 1
        theta_o = la0;
        ds_o = ds;
      end
    else
      % fprintf(1,'Lens angle range %.1f to %.1f\n', lens_angle(r,d,s,f,0), ...
      %   lens_angle(r,d,s,f,-r*d/(d^2+s^2)));
      x = [];
    end
  end
end
