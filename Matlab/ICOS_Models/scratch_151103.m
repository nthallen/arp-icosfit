%%
% Exploration of interleave options with R1!=R2, fixing L and determining
% realistic tolerances while adhering to the focusing and overlap
% criteria.
%
% This investigation qualifies as part of the sropt series, since it's
% working back to engineering a specific design from the constraints.
L = 50; % ICOS mirror spacing
B = 0.4; % Beam width, formerly known as 'W'
rd = .1; % detector radius (or half-width)
th = 15; % detector acceptance half-angle in degrees
C = 3000; % coherence length of laser
rdtanth = rd*tand(th);
mmin = ceil(C/(2*L));
plotnum = 2;
%%
clc;
for k=1:ceil(mmin/2)
  %%
  m = mmin;
  while(gcd(m,k) > 1)
    m = m + 1;
  end
  % R = L/(1-cos(k*pi/m)); % R1=R2 case
  %%
  rmin = B/(2*sin(pi/m));
  %%
  R1 = linspace(L,2000,2000)'; % arbitrary, probably want to improve
  R2 = L*(R1-L)./(R1*sin(k*pi/m)^2-L);
  z = [0; diff(sign(R2)) ~= 0];
  zi = (1:length(z))'+cumsum(z);
  zii = find(z)+(1:sum(z))'-1;
  R1(zi) = R1;
  R1(zii) = NaN;
  R2(zi) = R2;
  R2(zii) = NaN;
  dR1 = R1*.01;
  dR2 = R2*.01;
  R1d = R1*[.99 .99 1.01 1.01];
  R2d = R2*[.99 1.01 .99 1.01];
  % This section assumes we can adjust L to correct for R tolerance:
  Ld = ((R1d+R2d)+sqrt((R1d+R2d).^2-4.*R1d.*R2d*sin(k*pi/m)^2))/2;
  Lr = diff(minmax(Ld),1,2);
  plot(R1,Lr);
  xlabel('R1 cm');
  ylabel('L tolerance cm');
  ylim([0 10]);
  title(sprintf('L=%.0f Interleave %d/%d', L, k, m));
  %%
  % In this section, we hold L fixed and see if we can go far enough
  % out in r1 to make this work.
  Ll = L*ones(1,4);
  r1_2 = rdtanth*sqrt(((R2d-Ll).*R1d.^2.*Ll)./((R1d-Ll).*(R1d+R2d-Ll)));
%   r1_2(r1_2<0 | R1==L) = 0;
%   % r1_2(R1==L) = 0; % Skip the R1==L case
  r1max = sqrt(r1_2);
  dN = evaluate_interleave(m,k);
  phi = pi/m;
  phi_min = asin(B./(2*r1max));
  if any(phi_min >= phi)
    error('phi_min >= phi');
  end
  dphi = (phi-phi_min)./dN;
  philims = k*phi+dphi;
  sinphid = sqrt(Ll.*(R1d+R2d-Ll)./(R1d.*R2d));
%   if plotnum == 1
%     plot(R1,r1max,R1,rmin*ones(size(R1))); shg;
%     ylim([0 4]);
%     xlabel('R_1 cm');
%     ylabel('r cm');
%     legend('r_{max}', 'r_{min}','Location','SouthEast');
%     title(sprintf('Interleave %d/%d', k, m));
%   else
%     rvals = linspace(rmin,4,100);
%     [X,Y] = meshgrid(R1,rvals);
%     R2M = L*(X-L)./(X*sin(k*pi/m)^2-L);
%     s = Y.*sqrt((X-L).*(X+R2M-L)./(X.^2*L.*(R2M-L)));
%     sr = s.*Y;
%     %  mesh(X,Y,s); shg;
%     
%     wh1 = B/2;
%     Lh = wh1./s;
%     Lh(sr > rdtanth) = NaN;
%     
%     figure(1);
%     mesh(X,Y,Lh);
%     %   image(minmax(R1),minmax(rvals),Lh,'CDataMapping','scaled');
%     %   set(gca,'YDir','Normal');
%     xlabel('R_1 cm'); ylabel('r cm'); zlabel('Lh');
%     zlim([0 L]);
%     title(sprintf('L=%.0f Interleave %d/%d: min(Lh)=%.1f', L, k, m, ...
%       nanmin(nanmin(Lh))));
%     
%     figure(2);
%     v = any(~isnan(Lh));
%     plot(R1(v), R2(v));
%     xlabel('R_1 cm');
%     ylabel('R_2 cm');
%     ylim([-400 400]);
%     title(sprintf('L=%.0f Interleave %d/%d: min(Lh)=%.1f', L, k, m, ...
%       nanmin(nanmin(Lh))));
%   end
  pause;
end
