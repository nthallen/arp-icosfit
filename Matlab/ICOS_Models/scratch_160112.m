%%
% Review the skew * radius conjecture
load('IS_L50c.4.mat')
P = render_model(IS.res2(13));
P.ICOS_passes_per_injection = 30;
P.HR = 0;
% P.evaluate_endpoints = 2;
P.visible = 0;
PM = ICOS_Model6(P);
%%
[xyz, dxyz, d, s] = PM.M.extract_endpoints_skew(2);
as = abs(s);
r = sqrt(sum(xyz(:,[2,3]).^2,2));
msr = as(1).*r(1);
rel_sr = ((as.*r)-msr)/msr;
plot(rel_sr,'*-');
%%
for i=3:length(PM.M.Optic)
  [xyz, dxyz, d, s] = PM.M.extract_endpoints_skew(i);
  as = abs(s);
  r = sqrt(sum(xyz(:,[2,3]).^2,2));
  % msr = as(1).*r(1);
  rel_sr = ((as.*r)-msr)/msr;
  hold on;
  plot(rel_sr,'*-');
end
hold off;
