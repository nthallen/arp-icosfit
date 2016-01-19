%%
% Look at real sr evolution
load IS_L50c.4.mat
P = render_model(IS.res2(24), 'max_rays', 30000, 'visible', 0);
PM = ICOS_Model6(P);
%%
[oxyz, r, div, skew] = PM.M.extract_origin_skew(2);
sr = r.*abs(skew);
plot(sr,'.');
%%
sr0 = sr(1);
dsr = (sr-sr0)./sr0;
plot(dsr,'.');
%%
for opt_n = [3 4 5 6];
  [oxyz, r, div, skew] = PM.M.extract_origin_skew(opt_n);
  sr = r.*abs(skew);
  dsr = (sr-sr0)./sr0;
  plot(dsr,'.');
  title(sprintf('Optic %d', opt_n));
  pause;
end
