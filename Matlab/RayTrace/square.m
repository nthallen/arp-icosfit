classdef square < opt_surface
  properties
    l
    perimeter
    rad1
    rad2
    cos_theta
  end
  
  methods
    function sq = square(l, O, D, T, R, ni, ne, theta)
      % sq = square(l, O, D, T, R, ni, ne);
      % sq = square(l, O, D, T, R, ni, ne, theta);
      % l: linear size
      % O: origin vector
      % D: direction vector point out of external side
      % T: Transmittance
      % R: Reflectance
      % ni: Internal index of refraction
      % ne: External index of refraction
      % theta: Acceptance angle in degrees (optional, defaults to 90)
      %
      % Note: this square is missing rotation around the normal
      % vector. It must be extended when that flexibility is required.
      sq = sq@opt_surface(O, D, T, R, ni, ne);
      sq.l = l;
      
      if nargin < 8
        sq.cos_theta = 0;
      else
        sq.cos_theta = cosd(theta);
      end
      sq.rad1 = cross([0 0 1], D);
      nrad = sqrt(dot(sq.rad1,sq.rad1));
      if nrad < 0.1
        sq.rad1 = cross([0 1 0], D);
        nrad = sqrt(dot(sq.rad1,sq.rad1));
      end
      sq.rad1 = sq.rad1/nrad;
      sq.rad2 = cross(sq.rad1,D);
      
      l2 = l/2;
      col = ones(9,1);
      sq.perimeter = col*O + l2 * [0;1;1;1;0;-1;-1;-1;0]*sq.rad1 + ...
        l2 * [1; 1; 0; -1; -1; -1; 0; 1; 1]*sq.rad2;
    end
    
    function draw(sq)
      X = cell(3,1);
      for i=1:3
        X{i} = [ sq.perimeter(1:3,i)';
          [sq.perimeter(8,i) sq.O(i) sq.perimeter(4,i)];
          sq.perimeter([7 6 5],i)' ];
      end
      h = surfl(X{:});
      alpha(h,0.5);
      shading interp;
    end
    
    function [Pintercept,Vnormal] = intercept(sq, Rincident)
      if isempty(Rincident)
        Pintercept = [];
        Vnormal = [];
      elseif abs(dot(Rincident.D,sq.D)) < sq.cos_theta
        Pintercept = [];
        Vnormal = [];
      else
        M = [ sq.D 0
          1 0 0 -Rincident.D(1)
          0 1 0 -Rincident.D(2)
          0 0 1 -Rincident.D(3) ];
        A = [ dot(sq.O,sq.D); Rincident.O'];
        V = M\A; % V = [x y z t]'
        Pintercept = V(1:3)';
        Vnormal = sq.D;
        if V(4) < 0
          % does not intercept
          Pintercept = [];
        else
          dO = Pintercept - sq.O;
          if (abs(dot(dO,sq.rad1)) > sq.l/2) || ...
              (abs(dot(dO,sq.rad2)) > sq.l/2)
            Pintercept = [];
          end
        end
      end
    end
  end
end
