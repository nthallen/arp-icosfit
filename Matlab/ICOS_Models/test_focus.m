function test_focus(R, L, r)
% test_focus(R, L, r)
% Evaluates the possiblity of focusing the output of an ICOS cell
% with radius of curvature R, mirror spacing L and a circular spot
% pattern or radius r onto a detector of radius rd=0.l with an
% incident angle no greater than 15 degrees from normal.
rd = 0.1;
tan15 = tan(deg2rad(15));
rmax = sqrt(R*rd*tan15/sqrt(2*R/L-1));
fprintf(1,'Given R=%.0f and L=%.0f, you need r < %.1f to focus\n', R, L, rmax);
Lmin = (2*R*r^4)/(r^4+R^2*rd^2*tan15^2);
fprintf(1,'Given R=%.0f and r=%.1f, you need L > %.0f to focus\n', R, r, Lmin);

a = rd^2*tan15^2*L;
b = -2*r^4;
c = r^4*L;
disc = b^2-4*a*c;
if disc < 0
  if a > 0
    fprintf(1, 'Apparently no solutions for L=%.0f and r=%.1f\n', L, r);
  else
    fprintf(1, 'Apparently any R will work for L=%.0f and r=%.1f\n', L, r);
  end
else
  Rlims = sort((-b + [-1 1]*sqrt(disc))/(2*a));
  fprintf(1, 'Given L=%.0f and r=%.1f, you need ', L, r);
  if a > 0
    fprintf(1, 'R < %f or R > %f\n', Rlims);
  else
    fprintf(1, '%f <= R <= %f\n', Rlims);
  end
end

% n = 2.4361;
% a = (2*r^2)/(R*L) - r^2/R^2 + (r^2*n^2)/R^2 - 0.0718;
% b = -2*n*r^2/R;
% c = r^2;
% disc = b^2-4*a*c;
% if disc < 0
%   if a > 0
%     fprintf(1, 'Apparently no solutions (angle criteria)\n');
%   else
%     fprintf(1, 'Apparently any f will do (angle criteria)\n');
%   end
% else
%   f = sort((-b + [-1 1]*sqrt(disc))/(2*a));
%   if a > 0
%     fprintf(1, '%f <= f <= %f (angle criteria)\n', f);
%   else
%     fprintf(1, 'f < %f or f > %f (angle criteria)\n', f);
%   end
% end
% 
% a = 100*r^2*(2*R-L)-(2*R - L + n^2*L);
% b = 2*n*R*L;
% c = -R^2*L;
% % d = r/R;
% % s = d*sqrt((2*R-L)/L);
% % a = (100*r^2-1)*s^2*R^2;
% % b = 2*r^2*n*R;
% % c = -r^2*R^2;
% disc = b^2-4*a*c;
% if disc < 0
%   if a > 0
%     fprintf(1, 'Apparently no solutions (distance criteria)\n');
%   else
%     fprintf(1, 'Apparently any f will do (distance criteria)\n');
%   end
% else
%   f = sort((-b + [-1 1]*sqrt(disc))/(2*a));
%   if a > 0
%     fprintf(1, '%f <= f <= %f (distance criteria)\n', f);
%   else
%     fprintf(1, 'f < %f or f > %f (distance criteria)\n', f);
%   end
% end
% 
% if nargin > 3
%   % Test calculation for given values of ftest
%   d = r/R;
%   fprintf(1, '\nd = %f\n', d);
%   s = d*sqrt((2*R-L)/L);
%   fprintf(1, 's = %f\n\n', s);
%   for i = 1:length(ftest)
%     f = ftest(i);
%     fprintf(1, 'For f = %f:\n', f);
%     d2 = r*(n/R - 1/f);
%     fprintf(1, '  d2 = %f\n', d2);
%     xmin = -(r*d2)/(d2^2+s^2);
%     fprintf(1, '  xmin = %f\n', xmin);
%     rx = r*s/sqrt(s^2+d2^2);
%     fprintf(1, '  rx = %f\n', rx);
%     sx = sqrt(s^2+d2^2);
%     fprintf(1, '  sx = %f\n', sx);
%     theta = rad2deg(atan(sx));
%     fprintf(1, '  theta = %f\n', theta);
%   end
% end
