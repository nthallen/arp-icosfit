classdef detector < optic
  properties
    l % side length
  end
  
  methods
    function det = detector(l, O, D, vis, theta)
      % det = detector(l, O, D, vis, theta);
      % Models a square detector element
      % l: linear size
      % O: origin vector
      % D: direction vector
      % vis: visability boolean
      % theta: [Optional] Acceptance angle in degrees
      if nargin < 5
        theta = 0;
      end
      det = det@optic('det', O, D, 1, 1, vis);
      det.Surface{1} = square(l, O, D, 0, 0, 1, 1, theta);
      det.Surface{2} = square(l, O-0.1*D, -D, 0, 0, 1, 1);
    end
  end
end
