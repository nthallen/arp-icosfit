classdef spherical_segment < opt_surface
  properties
    r % radius
    RC % radius of curvature
    OC % origin of curvature
    perimeter
  end
  
  properties (Access = private)
    rad1
    rad2
    mincosphi
  end
  
  methods
    function cap = spherical_segment(r, RC, O, D, T, R, ni, ne)
      % cap = spherical_surface(r, RC, O, D, T, R, ni, ne)
      % r: radius
      % RC: radius of curvature (- = convex, + = concave)
      % O: origin
      % D: unit normal pointing out
      % T: Transmittance
      % R: Reflectance
      % ni: internal index of refraction
      % ne: external index of refraction
      cap = cap@opt_surface(O, D, T, R, ni, ne);
      cap.r = r;
      cap.RC = RC;
      cap.OC = O + RC*D;
      res = 50;
      cap.rad1 = cross([0 0 1], D);
      nrad = sqrt(dot(cap.rad1,cap.rad1));
      if nrad < 0.1
        cap.rad1 = cross([0 1 0], D);
        nrad = sqrt(dot(cap.rad1,cap.rad1));
      end
      cap.rad1 = cap.rad1/nrad;
      cap.rad2 = cross(cap.rad1,D);
      res = 50;
      th = linspace(0,2*pi,res)';
      col = ones(length(th),1);
      phimax = asin(cap.r/cap.RC);
      cap.mincosphi = cos(phimax);
      cap.perimeter = ...
        col*cap.OC -cap.RC*cap.mincosphi*col*cap.D + r*sin(th)*cap.rad1 + r*cos(th)*cap.rad2;
    end
    
    function draw(cap)
      res = 50;
      th = linspace(0,2*pi,res)';
      phimax = asin(cap.r/cap.RC);
      ph = linspace(0,phimax,10)';
      [phi,theta] = meshgrid(ph,th);
      u = cap.RC*cos(phi);
      v = cap.RC*sin(phi).*sin(theta);
      w = cap.RC*sin(phi).*cos(theta);
      X = cap.OC(1) - u*cap.D(1) + v*cap.rad1(1) + w*cap.rad2(1);
      Y = cap.OC(2) - u*cap.D(2) + v*cap.rad1(2) + w*cap.rad2(2);
      Z = cap.OC(3) - u*cap.D(3) + v*cap.rad1(3) + w*cap.rad2(3);
      h = surfl(X,Y,Z);
      alpha(h,0.5);
      shading interp;
    end
    
    function [Pintercept,Vnormal] = intercept(cap, Rincident)
      if isempty(Rincident)
        Pintercept = [];
        Vnormal = [];
      else
        dO = Rincident.O-cap.OC;
        a = sum(Rincident.D.^2);
        b = 2*sum(dO.*Rincident.D);
        c = sum(dO.^2) - cap.RC^2;
        disc = b^2 - 4*a*c;
        if disc <= 0 % we'll call that a miss
          t = [];
        else
          t = (-b + [1,-1]*sqrt(disc))/(2*a);
          t = t(t > 0);
          if length(t) > 1
            t = min(t);
          end
        end
        if isempty(t)
          Pintercept = [];
          Vnormal = [];
        else
          Pintercept = Rincident.O + t*Rincident.D;
          cosphi = dot((Pintercept-cap.OC)/cap.RC,(cap.O-cap.OC)/cap.RC);
          if cosphi < cap.mincosphi
            Pintercept = [];
            Vnormal = [];
          else
            Vnormal = cap.OC - Pintercept;
            if cap.RC < 0
              Vnormal = -Vnormal;
            end
            Vnormal = Vnormal/sqrt(sum(Vnormal.^2));
          end
        end
      end
    end
  end
end
