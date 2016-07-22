classdef HRmirrorP < optic
  % HRmirrorP: Pressure distorted plano-concave or -convex mirror
  properties
    % These properties could be simply passed on to the surfaces
    r % mirror radius
    RC % radius of curvature
    CT % mirror center thickness
    T % transmittance
    R % reflectance
    dP % differential pressure (default: 0 Torr)
    E % Young's Modulus (default: 10.1e6 for ZnSe)
    pr % Poisson's ratio (default: 0.25 for ZnSe)
  end
  
  methods
    function HR = HRmirrorP(nm, r_in, RC_in, CT, T_in, R_in, O_in, D_in, ...
        ni_in, ne_in, vis, varargin)
      % HR = HRmirror(name, r, RoC, w, T, R, O, D, ni, ne, visible, ...
      %   [dP [, E[ , pr]]]);
      % name: name
      % r: mirror radius, cm
      % RoC: radius of curvature, cm
      % w: center thickness, cm
      % T: transmitance of the curved surface
      % R: reflectance of the curved surface
      % O: position of the center of the curved surface
      % D: direction outward from O
      % dP: differential pressure in direction of D in torr
      % E: Young's Modulus
      % pr: Poisson's ratio
      % ni: index of refraction of the mirror
      % ne: index of refraction of the air
      % vis: boolean indicating whether or not to draw the mirror
      HR = HR@optic(nm, O_in, D_in, ni_in, ne_in, vis);
      HR.r = r_in;
      HR.RC = RC_in;
      HR.CT = CT;
      HR.T = T_in;
      HR.R = R_in;
      HR.dP = 0;
      HR.E = 10.1e6;
      HR.pr = 0.28;
      Ocirc = HR.O - CT*HR.D;
      Dcirc = -HR.D;
      if length(varargin) >= 1
        HR.dP = varargin{1};
        if length(varargin) >= 2
          HR.E = varargin{2};
          if length(varargin) >= 3
            HR.pr = 0.28;
          end
        end
      end
      k = 3*(HR.dP*14.7/760)*(1-HR.pr^2)/(16*HR.E*(CT^3)*2.54);
      HR.Surface{1} = spherical_segmentP(HR.r, HR.RC, HR.O, HR.D, ...
        T_in, R_in, ni_in, ne_in, k);
      HR.Surface{2} = spherical_segmentP(HR.r, 0, Ocirc, Dcirc, 1, 0, ...
        ni_in, ne_in, -k);
    end
  end
end
