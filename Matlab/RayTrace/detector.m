classdef detector < optic
  properties
    l % side length
  end
  
  methods
    function det = detector(l, O, D, vis)
      det = det@optic('det', O, D, 1, 1, vis);
      det.Surface{1} = square(l, O, D, 0, 0, 1, 1);
      det.Surface{2} = square(l, O-0.1*D, -D, 0, 0, 1, 1);
    end
  end
end