classdef square_window < optic
  properties
    l % side length
  end
  
  methods
    function win = square_window(l, th, O, D, vis)
      % win = square_window(l, th, O, D, vis);
      % l: linear size
      % th: thickness
      % O: origin vector
      % D: direction vector
      % vis: visability boolean
      % Models a square window element
      win = win@optic('sqwin', O, D, 1, 1, vis);
      win.Surface{1} = square(l, O, D, 1, 0, 1, 1);
      win.Surface{2} = square(l, O-th*D, -D, 1, 0, 1, 1);
    end
  end
end
