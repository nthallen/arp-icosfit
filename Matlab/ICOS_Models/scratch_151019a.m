%%
% Exploration of interleave options with R1=R2, fixing L and determining
% realistic tolerances while adhering to the focusing and overlap
% criteria.
L = 50;
W = 0.4;
rd = .1;
th = 15;
C = 3000;
rdtanth = rd*tand(th);
mmin = ceil(C/(2*L));
%%
clc;
for k=1:ceil(mmin/2)
  %%
  m = mmin;
  while(gcd(m,k) > 1)
    m = m + 1;
  end
  R = L/(1-cos(k*pi/m));
  %%
  rmin = W/(2*sin(pi/m));
  %%
  s = (rmin/R)*sqrt((2*R/L)-1);
  sr = s*rmin;
  if R < 2000 && sr < rdtanth
    rmax = sqrt(rdtanth*R/sqrt((2*R/L)-1));
    fprintf(1, 'Interleave: %d/%d R: %.2f sr/rdtanth: %.3f  r: [%f,%f]\n', ...
      m, k, R, sr/rdtanth, rmin, rmax);
    dN = evaluate_interleave(m,k);
    phi = pi/m;
    phi_min = asin(W/(2*rmax));
    if phi_min >= phi
      error('phi_min >= phi');
    end
    dphi = (phi-phi_min)./dN;
    philims = k*phi+dphi;
    Rlims = L./(1-cos(philims));
    sag = Rlims - sqrt(Rlims.^2-rmax^2);
    dsag = abs(diff(sag))*10; % convert to mm
    fprintf(1,'  R in [%.1f, %.1f] cm, dsag = %.3f mm\n', Rlims(2), Rlims(1), dsag);
    Llims = R.*(1-cos(philims));
    fprintf(1,'  L in [%.3f, %.3f] cm, dL = %.3f mm\n', Llims(2), Llims(1), ...
      abs(diff(Llims))*10);
    % Now suppose we can buy mirrors with a +/- 10 um sag tolerance
    % What range of R would that give us, and what range of L would
    % we need to accomodate that and still meet our angle target
    rmirror = 1.5*2.54; % 3" mirror (conservatively)
    sags = R - sqrt(R^2-rmirror^2) + [-1 1]*5e-4;
    L1 = R*(1-cos(k*phi));
    Rs = (sags.^2 + rmirror^2)./(2*sags);
    Ls = Rs*(1-cos(k*phi));
    fprintf(1,'  R mfr [%.1f, %.1f]  L range: [%.1f, %.1f] kphi = %.1f\n', ...
      Rs(1), Rs(2), Ls(1), Ls(2), k*phi);
  end
end
%%
% Look at tolerances on short R
% W/(2*rmax) <= w/r <= sin(pi/n)
rmax = 3*2.54/2; % 3" diameter max
n = ceil(C/(2*L));
phi = pi/n;
t = (W/(2*rmax))^2;
R0 = L*(1-sqrt(1-t))/t;
R1 = L/(1+cos(phi));
sag0 = R0-sqrt(R0^2-rmax^2);
sag1 = R1-sqrt(R1^2-rmax^2);
dsag = sag0-sag1;
