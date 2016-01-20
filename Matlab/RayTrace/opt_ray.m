classdef opt_ray
  properties
    O % origin (row)
    D % direction (row)
    E % endpoint
    P % Power
    pass % pass number
    Inside % boolean
  end
  methods
    function R = opt_ray(O, D)
      % R = opt_ray(O, D);
      R.O = O;
      R.D = D/sqrt(sum(D.^2));
      R.E = O+D;
      R.P = 1;
      R.pass = 1;
      R.Inside = 1;
    end
    
    function [h_out, th_out] = draw(R, varargin)
      x = [R.O(1) R.E(1)];
      y = [R.O(2) R.E(2)];
      z = [R.O(3) R.E(3)];
      h = plot3(x, y, z, 'r');
      if ~isempty(varargin)
        th = mtext(mean(x), mean(y), mean(z), varargin{:});
      end
      if nargout > 0
        h_out = h;
        if nargout > 1
          th_out = th;
        end
      end
    end
    
    function [srO,srE] = sr(R)
      srO = (-R.O(3)*R.D(2)+R.O(2)*R.D(3))/R.D(1);
      if nargout > 1
        srE = (-R.E(3)*R.D(2)+R.E(2)*R.D(3))/R.D(1);
      end
    end
    
    function [DL,L] = len(R)
      DL = sqrt(sum(R.D.^2));
      if nargout > 1
        L = sqrt(sum((R.E-R.O).^2));
      end
    end
  end
end
