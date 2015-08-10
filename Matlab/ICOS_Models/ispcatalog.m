function Cat = ispcatalog
% Cat = ispcatalog;
% Builds a struct array of plano-convex ZiSe blanks that
% we could use for high-reflectivity mirrors
Cat = [
  ispoptic('Custom75', 3, -750)
  ispoptic('ZC-PX-12-25', 0.5, 35.97)
  ispoptic('ZC-PX-12-50', 0.5, 71.46)
  ispoptic('ZC-PX-25-50', 1, 71.46)
  ispoptic('ZC-PX-38-50', 1.5, 71.46)
  ispoptic('ZC-PX-25-63', 1, 89.14)
  ispoptic('ZC-PX-38-63', 1.5, 89.14)
  ispoptic('ZC-PX-25-76', 1, 106.15)
  ispoptic('ZC-PX-50-76', 2, 106.15)
  ispoptic('ZC-PX-38-76', 1.5, 106.15)
  ispoptic('ZC-PX-25-100', 1, 139.96)
  ispoptic('ZC-PX-38-100', 1.5, 139.96)
  ispoptic('ZC-PX-50-100', 2, 139.96)
  ispoptic('ZC-PX-25-127', 1, 180.32)
  ispoptic('ZC-PX-38-127', 1.5, 180.32)
  ispoptic('ZC-PX-50-127', 2, 180.32)
  ispoptic('ZC-PX-25-150', 1, 209.91)
  ispoptic('ZC-PX-38-150', 1.5, 209.91)
  ispoptic('ZC-PX-50-150', 2, 209.91)
  ispoptic('ZC-PX-25-200', 1, 281.20)
  ispoptic('ZC-PX-38-200', 1.5, 281.20)
  ispoptic('ZC-PX-50-200', 2, 281.20)
  ispoptic('ZC-PX-25-254', 1, 356.51)
  ispoptic('ZC-PX-38-254', 1.5, 356.51)
  ispoptic('ZC-PX-50-254', 2, 356.51)
  ispoptic('ZC-PX-50-300', 2, 424.67)
  ispoptic('ZC-PX-25-500', 1, 698.17)
  ispoptic('ZC-PX-38-500', 1.5, 698.17)
  ispoptic('ZC-PX-50-500', 2, 698.17)
  ispoptic('ZC-PX-38-1000', 1.5, 1404.00)
  ispoptic('ZC-PX-50-1000', 2, 1404.00)
  ispoptic('ZC-PX-38-2000', 1.5, 2808.00)
  ];

function elt = ispoptic(catnum, dia, R)
% input radius of curvature is in mm. Output is in cm.
elt = struct('cat', catnum, 'dia_in', dia, 'R_cm', -R/10);
