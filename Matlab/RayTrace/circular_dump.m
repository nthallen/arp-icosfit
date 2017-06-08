classdef circular_dump < optic
  properties
    l % side length
  end
  
  methods
    function win = circular_dump(r, th, O, D, vis)
      % win = circular_dump(r, th, O, D, vis);
      % r: radius
      % th: thickness
      % O: origin vector
      % D: direction vector
      % vis: visability boolean
      % Models a circular beam dump element
      win = win@optic('circwin', O, D, 1, 1, vis);
      win.Surface{1} = circular_surface(r, O, D, 0, 0, 1, 1);
      win.Surface{2} = circular_surface(r, O-th*D, -D, 0, 0, 1, 1);
    end
  end
end
