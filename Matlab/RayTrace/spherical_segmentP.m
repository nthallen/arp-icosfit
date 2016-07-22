classdef spherical_segmentP < opt_surface
  properties
    r % radius
    RC % radius of curvature
    perimeter
    k % pressure coefficient
  end
  
  properties (Access = private)
    rad1
    rad2
    mincosphi
  end
  
  methods
    function cap = spherical_segmentP(r, RC, O, D, T, R, ni, ne, k)
      % cap = spherical_segmentP(r, RC, O, D, T, R, ni, ne, k)
      % r: radius
      % RC: radius of curvature (- = convex, + = concave)
      % O: origin
      % D: unit normal pointing out
      % T: Transmittance
      % R: Reflectance
      % ni: internal index of refraction
      % ne: external index of refraction
      % k: pressure perturbation coefficient
      cap = cap@opt_surface(O, D, T, R, ni, ne);
      cap.r = r;
      cap.k = k;
      cap.RC = RC;
      cap.O = O;
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
      if cap.RC ~= 0
        phimax = asin(cap.r/cap.RC);
        cap.mincosphi = cos(phimax);
        cap.perimeter = ...
          col*cap.O +cap.RC*(1-cap.mincosphi*col)*cap.D + ...
            r*sin(th)*cap.rad1 + r*cos(th)*cap.rad2;
      else
        cap.perimeter = col*cap.O + r*sin(th)*cap.rad1 + r*cos(th)*cap.rad2;
      end
    end
    
    function draw(cap)
      res = 50;
      th = linspace(0,2*pi,res)';
      rv = linspace(0,cap.r,10)';
      %phimax = asin(cap.r/cap.RC);
      %ph = linspace(0,phimax,10)';
      [rm,theta] = meshgrid(rv,th);
      % u = cap.RC*cos(phi);
      if cap.RC == 0
        u = cap.k*(cap.r^2-rm.^2).^2;
      else
        u = cap.RC - sqrt(cap.RC^2-rm.^2) + cap.k*(cap.r^2-rm.^2).^2;
      end
      v = rm.*sin(theta);
      w = rm.*cos(theta);
      X = cap.O(1) - u*cap.D(1) + v*cap.rad1(1) + w*cap.rad2(1);
      Y = cap.O(2) - u*cap.D(2) + v*cap.rad1(2) + w*cap.rad2(2);
      Z = cap.O(3) - u*cap.D(3) + v*cap.rad1(3) + w*cap.rad2(3);
      h = surfl(X,Y,Z);
      alpha(h,0.5);
      shading interp;
    end
    
    function [Pintercept,Vnormal] = intercept(cap, Rincident)
      if isempty(Rincident)
        Pintercept = [];
        Vnormal = [];
      else
        % X(t) is position relative to surface origin at distance t along
        % the ray.
        X = @(t) Rincident.O+t*Rincident.D-cap.O;
        % dxR(t) is the component of the relative position normal to the
        % mirror's center surface.
        dxR = @(t) sum(X(t).*cap.D);
        % rv(t) is the radial vector from the optical axis to X(t)
        rv = @(t) X(t)-dxR(t)*cap.D;
        % r2(t) is the square of the length of the radius
        r2 = @(t) sum(rv(t).^2);
        if cap.RC == 0
          dxr = @(t) cap.k*(cap.r^2-r2(t))^2;
        else
          dxr = @(t) cap.RC - sqrt(cap.RC^2-r2(t)) + cap.k*(cap.r^2-r2(t))^2;
        end
        ddxr = @(t) dxR(t)-dxr(t);
        t = fzero(ddxr,0);
        r2t = r2(t);
        if r2t > cap.r^2
          t = [];
        end
        if isempty(t)
          Pintercept = [];
          Vnormal = [];
        else
          Pintercept = Rincident.O + t*Rincident.D;
          rvt = rv(t);
          ri = sqrt(r2t);
          if cap.RC == 0
            drdx = 4*ri*cap.k*(cap.r^2-r2t);
          else
            drdx = - ri/sqrt(cap.RC^2-r2t) + 4*ri*cap.k*(cap.r^2-r2t);
          end
          if ri == 0
            Vnormal = cap.D;
          else
            Vnormal = rvt*drdx/ri + cap.D;
            Vnormal = Vnormal/sqrt(sum(Vnormal.^2));
          end
        end
      end
    end
  end
end
