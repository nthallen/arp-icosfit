classdef opt_surface
  % opt_surface
  % Optical surface abstract base class
  properties
    O
    D % direction pointing to the external side
    T % transmittance
    R % reflectance
    n_int
    n_ext
    visible
    emission_threshold
  end
  
  properties (Abstract)
    perimeter
  end
  
  methods
    function surf = opt_surface(O, D, T, R, n_int, n_ext)
      surf.O = O;
      surf.D = D;
      surf.T = T;
      surf.R = R;
      surf.n_int = n_int;
      surf.n_ext = n_ext;
      surf.visible = true;
      surf.emission_threshold = 0;
    end
    
    function [Rincident,Rref,Rtrans] = propagate(surf, Rincident)
      [Pintercept,Vnormal] = surf.intercept(Rincident);
      % Pintercept is Rincident extended to the intersection point
      if isempty(Pintercept)
        if ~isempty(Rincident) && surf.visible
          Rincident.draw;
        end
        Rref = [];
        Rtrans = [];
      else
        Rincident.E = Pintercept;
        if surf.visible
          Rincident.draw;
        end
        N = dot(Rincident.D,Vnormal);
        if surf.R == 0
          Rref = [];
        else
          Rref = opt_ray(Rincident.E, Rincident.D - 2*N*Vnormal);
          Rref.P = Rref.P * surf.R;
          if Rref.P < surf.emission_threshold
            Rref = [];
          end
        end
        if surf.T == 0
          Rtrans = [];
        else
          % handle refraction
          sininc = Rincident.D - N*Vnormal;
          if N > 0
            sintrans = sininc * surf.n_int / surf.n_ext;
          else
            sintrans = sininc * surf.n_ext / surf.n_int;
          end
          Rtrans = opt_ray(Pintercept, N*Vnormal+sintrans);
          Rtrans.P = Rincident.P * surf.T;
          if Rtrans.P < surf.emission_threshold
            Rtrans = [];
          end
        end
      end
    end
  end
  
  methods (Abstract)
    [Pintercept,Vnormal] = intercept(surf, Rincident)
    draw(surf)
  end
end