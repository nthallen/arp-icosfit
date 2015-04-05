classdef Herriott_Mirror < HRmirror
  properties
    aperature_point
    aperature_radius
  end
  
  methods
    function HM = Herriott_Mirror(nm, r_in, RC_in, CT_in, R_in, Ap, Ar, O_in, D_in, vis)
      % HM = Herriott_Mirror(name, r, RoC, CT, R, O, D, ni, ne, visible);
      % name: name
      % r: mirror radius
      % RoC: radius of curvature
      % CT: center thickness of mirror
      % R: Mirror reflectivity
      % Ar: Aperature radius (beam diameter)
      % O: position of the center of the curved surface
      % D: direction from O toward center of curvature
      % vis: boolean indicating whether or not to draw the mirror
      HM = HM@HRmirror(nm, r_in, RC_in, CT_in, 0, R_in, O_in, D_in, 1, 1, vis);
      HM.aperature_point = Ap;
      HM.aperature_radius = Ar;
    end
    
    function [Rincident, Rreflect, Rinternal, Rtransmit] = propagate(HM, Rincident)
      [Rincident, Rreflect, Rinternal, Rtransmit] = HM.propagate@HRmirror(Rincident);
      if isempty(HM.aperature_point)
        HM.aperature_point = Rincident.E;
        Rtransmit = Rincident;
        Rtransmit.O = Rtransmit.E;
        Rreflect = [];
        Rinternal = [];
      else
        da = Rincident.E - HM.aperature_point;
        lda = dot(da, HM.D);
        rda = da - lda*HM.D;
        r = sqrt(sum(rda.^2));
        if r < HM.aperature_radius
          d = dot(Rincident.E-Rincident.O,Rincident.D);
          Rincident.E = Rincident.O + (d+1)*Rincident.D;
          if HM.visible
            Rincident.draw;
          end
          Rreflect = [];
          Rtransmit = [];
        end
      end
    end
  end
end
