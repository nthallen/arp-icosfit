%%
P = ICOS_Model2.props;
P.RC = 150;
P.vis = true;
PM = ICOS_Model2(P);
%%
dz = .01:.005:.03;
rskew = 0*dz;
for i = 1:length(dz);
  P.dz = dz(i);
  P.injection_scale = 0.5;
  P.vis = false;
  PM = ICOS_Model2(P);
  xyz = PM.M.extract_endpoints(2);
  rmax = max(sqrt(sum(xyz(:,[2 3]).^2,2)));
  P.injection_scale = P.injection_scale * 3.5 / rmax;
  P.vis = true;
  PM = ICOS_Model2(P);
  [xyz,dxyz,div,skew] = PM.M.extract_endpoints_skew(2);
  r = sqrt(sum(xyz(:,[2 3]).^2,2));
  rskew(i) = mean(r.*skew);
end
figure;
plot(dz, rskew,'*');
xlabel('dz');
ylabel('r * skew');
