classdef Herriott_Mirror < HRmirror
  properties
    aperature_point
    aperature_radius
  end
  
  methods
    function HM = Herriott_Mirror(nm, r_in, RC_in, CT_in, R_in, Ap, Ar, O_in, D_in, vis)
      % HM = Herriott_Mirror(name, r, RoC, CT, R, Ap, Ar, O, D, ni, ne, visible);
      % name: name
      % r: mirror radius
      % RoC: radius of curvature
      % CT: center thickness of mirror
      % R: Mirror reflectivity
      % Ap: Aperature position
      % Ar: Aperature radius (beam diameter)
      % O: position of the center of the curved surface
      % D: direction from O toward center of curvature
      % vis: boolean indicating whether or not to draw the mirror
      HM = HM@HRmirror(nm, r_in, RC_in, CT_in, 0, R_in, O_in, D_in, 1, 1, vis);
      HM.aperature_point = Ap;
      HM.aperature_radius = Ar;
      % Surface{1} is spherical
      % Surface{2} is planar
      % Append Herriott aperature to the perimeter for Surface{1}
      res = 20;
      th = linspace(0,2*pi,res)';
      col = ones(length(th),1);
      Ap_peri = [col*Ap(1), Ap(2)+Ar*cos(th), Ap(3)+Ar*sin(th)];
      HM.Surface{1}.perimeter = [
        HM.Surface{1}.perimeter
        NaN NaN NaN
        Ap_peri
        ];
    end
    
    function [Rincident, Rreflect, Rinternal, Rtransmit] = propagate(HM, Rincident)
      [Rincident, Rreflect, Rinternal, Rtransmit] = HM.propagate@HRmirror(Rincident);
      if ~isempty(Rincident)
        if isempty(HM.aperature_point)
          HM.aperature_point = Rincident.E;
          Rtransmit = Rincident;
          Rtransmit.O = Rtransmit.E;
          Rreflect = [];
          Rinternal = [];
        else
          try
            da = Rincident.E - HM.aperature_point;
          catch
            error('Something broke');
          end
          lda = dot(da, HM.D);
          rda = da - lda*HM.D;
          r = sqrt(sum(rda.^2));
          if r < HM.aperature_radius
            d = dot(Rincident.E-Rincident.O,Rincident.D);
            Rincident.E = Rincident.O + (d+1)*Rincident.D;
            Rincident.Inside = -1; % Not inside, but not "outside"
            if HM.visible
              Rincident.draw;
            end
            Rreflect = [];
            Rtransmit = [];
          end
        end
        if HM.alternate && Rincident.Inside < 1
          Rincident = [];
        end
      end
    end
  end
end
