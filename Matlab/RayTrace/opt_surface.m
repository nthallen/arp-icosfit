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
    % emission_threshold
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
      % surf.emission_threshold = 0;
    end
    
    function [Rincident,Rref,Rtrans] = propagate(surf, Rincident)
      [Pintercept,Vnormal] = surf.intercept(Rincident);
      % Pintercept is Rincident extended to the intersection point
      if isempty(Pintercept)
        if ~isempty(Rincident)
          Pintercept = surf.intercept_plane(Rincident);
          if ~isempty(Pintercept)
            Rincident.E = Pintercept;
            Rincident.Inside = 0;
          end
          if surf.visible
            Rincident.draw;
          end
        end
        Rref = [];
        Rtrans = [];
      else
        Rincident.E = Pintercept;
        if surf.visible
          Rincident.draw;
        end
        N = dot(Rincident.D,Vnormal);
        if surf.R == 0 % || Rincident.P < surf.emission_threshold
          Rref = [];
        else
          Rref = opt_ray(Rincident.E, Rincident.D - 2*N*Vnormal);
          Rref.P = Rincident.P * surf.R;
          Rref.pass = Rincident.pass;
%           if Rref.P < surf.emission_threshold
%             Rref = [];
%           end
        end
        if surf.T == 0 % || Rincident.P < surf.emission_threshold
          Rtrans = [];
        else
          % handle refraction
          sininc = Rincident.D - N*Vnormal;
          if N > 0
            sintrans = sininc * surf.n_int / surf.n_ext;
          else
            sintrans = sininc * surf.n_ext / surf.n_int;
          end
          % Rtrans = opt_ray(Pintercept, N*Vnormal+sintrans);
          Rtrans = opt_ray(Pintercept, sign(N)*sqrt(1-sum(sintrans.^2))*Vnormal+sintrans);
          Rtrans.P = Rincident.P * surf.T;
          Rtrans.pass = Rincident.pass;
%           if Rtrans.P < surf.emission_threshold
%             Rtrans = [];
%           end
        end
      end
    end
    
    function Pintercept = intercept_plane(surf, Rincident)
      if isempty(Rincident)
        Pintercept = [];
      else
        M = [ surf.D 0
          1 0 0 -Rincident.D(1)
          0 1 0 -Rincident.D(2)
          0 0 1 -Rincident.D(3) ];
        A = [ dot(surf.O,surf.D); Rincident.O'];
        V = M\A; % V = [x y z t]'
        if V(4) < 0
          Pintercept = [];
        else
          Pintercept = V(1:3)';
        end
      end
    end
  end
  
  methods (Abstract)
    [Pintercept,Vnormal] = intercept(surf, Rincident)
    draw(surf)
  end
end
