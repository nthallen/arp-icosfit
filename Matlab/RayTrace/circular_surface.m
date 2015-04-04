classdef circular_surface < opt_surface
  properties
    r % radius
    perimeter
  end
  
  methods
    function circ = circular_surface(r, O, D, T, R, ni, ne)
      % circ = circular_surface(r, O, D, T, R, ni, ne)
      % r: radius
      % O: origin
      % D: unit normal pointing out
      % T: Transmittance
      % R: Reflectance
      % ni: internal index of refraction
      % ne: external index of refraction
      circ = circ@opt_surface(O, D, T, R, ni, ne);
      circ.r = r;
      res = 50;
      rad1 = cross([0 0 1], D);
      nrad = sqrt(dot(rad1,rad1));
      if nrad < 0.1
        rad1 = cross([0 1 0], D);
        nrad = sqrt(dot(rad1,rad1));
      end
      rad1 = rad1/nrad;
      rad2 = cross(rad1,D);
      th = linspace(0,2*pi,res)';
      col = ones(res,1);
      circ.perimeter = col*O + r*sin(th)*rad1 + r*cos(th)*rad2;
    end
    
    function draw(circ)
      col = ones(size(circ.perimeter,1),1);
      X = [ circ.perimeter(:,1) (circ.perimeter(:,1)+col*circ.O(1))/2 col*circ.O(1) ];
      Y = [ circ.perimeter(:,2) (circ.perimeter(:,2)+col*circ.O(2))/2 col*circ.O(2) ];
      Z = [ circ.perimeter(:,3) (circ.perimeter(:,3)+col*circ.O(3))/2 col*circ.O(3) ];
      h = surfl(X,Y,Z);
      alpha(h,0.5);
      shading interp;
    end
    
    function [Pintercept,Vnormal] = intercept(circ, Rincident)
      Pintercept = circ.intercept_plane(Rincident);
      if isempty(Pintercept)
        Vnormal = [];
      elseif sum((Pintercept-circ.O).^2) > circ.r^2
        Pintercept = [];
        Vnormal = [];
      else
        Vnormal = circ.D;
      end
    end
  end
end
