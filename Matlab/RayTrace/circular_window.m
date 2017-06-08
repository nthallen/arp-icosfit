classdef circular_window < optic
  properties
    l % side length
  end
  
  methods
    function win = circular_window(r, th, O, D, vis)
      % win = circular_window(r, th, O, D, vis);
      % r: radius
      % th: thickness
      % O: origin vector
      % D: direction vector
      % vis: visability boolean
      % Models a circular window element
      win = win@optic('circwin', O, D, 1, 1, vis);
      win.Surface{1} = circular_surface(r, O, D, 1, 0, 1, 1);
      win.Surface{2} = circular_surface(r, O-th*D, -D, 1, 0, 1, 1);
    end
  end
end
