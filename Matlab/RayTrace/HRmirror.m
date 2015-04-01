classdef HRmirror < optic
  properties
    % These properties could be simply passed on to the surfaces
    r % mirror radius
    RC % radius of curvature
    w % mirror width at center
    T % transmittance
    R % reflectance
  end
  
  methods
    function HR = HRmirror(nm, r_in, RC_in, w, T_in, R_in, O_in, D_in, ni_in, ne_in, vis)
      % HR = HRmirror(name, r, RoC, w, T, R, O, D, ni, ne, visible);
      % name: name
      % r: mirror radius
      % RoC: radius of curvature
      % w: width at mirror center
      % T: transmitance of the curved surface
      % R: reflectance of the curved surface
      % O: position of the center of the curved surface
      % D: direction from O toward center of curvature
      % ni: index of refraction of the mirror
      % ne: index of refraction of the air
      % vis: boolean indicating whether or not to draw the mirror
      HR = HR@optic(nm, O_in, D_in, ni_in, ne_in, vis);
      HR.r = r_in;
      HR.RC = RC_in;
      HR.T = T_in;
      HR.R = R_in;
      Ocirc = HR.O - w*HR.D;
      Dcirc = -HR.D;
      HR.Surface{1} = spherical_segment(HR.r, HR.RC, HR.O, HR.D, T_in, 1-T_in, ni_in, ne_in);
      HR.Surface{2} = circular_surface(HR.r, Ocirc, Dcirc, 1, 0, ni_in, ne_in);
      %HR.Surface{1} = circular_surface(HR.r, HR.O, HR.D, 1, 0, ni_in, ne_in);
    end
  end
end