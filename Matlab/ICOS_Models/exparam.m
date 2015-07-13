function R = exparam(P)
% R = exparam(P)
% P is a struct defining R1, r1, r2 and L (and Rw1).
% or P defines R1, R2, RR1, r1 and Rw1
% One of these may be a vector.
n = 2.4361;

R.RR1 = [];
R.Rr1 = [];
R.Rh1 = [];
R.Rw1 = [];
R.Rd1 = [];
R.Rs1 = [];
R.RL = [];
R.RR2 = [];
R.Rr2 = [];
R.Rh2 = [];
R.Rw2 = [];
R.Rd2 = [];
R.Rs2 = [];

R.R1 = [];
R.r1 = [];
R.h1 = [];
R.w1 = [];
R.d1 = [];
R.s1 = [];
R.L = [];
R.R2 = [];
R.r2 = [];
R.h2 = [];
R.w2 = [];
R.d2 = [];
R.s2 = [];

R.r15 = [];

flds = fields(P);
for i=1:length(flds)
  fld = flds{i};
  R.(fld) = P.(fld);
end

if isfield(P,'RR1') % assume RR1, R1, R2, r1 and Rw1
  R.d1 = R.r1./R.R1;
  R.Rd2 = -R.d1*n;
  R.RR2 = - R.R1/n;
  R.Rr2 = R.r1;
  V = [
    -R.Rd2.^2
    R.Rd2.^2.*R.RR1+2*R.Rr2.*R.Rd2
    -R.Rr2.*R.Rd2.*R.RR1 - R.Rr2.^2 - R.Rw1.^2
    R.Rw1.^2.*R.RR1
    ];
  rts = roots(V);
  rts = rts(imag(rts)==0);
  rts = rts(rts > 0);
  R.RL = rts'; % make it a row.
  R.Rh1 = R.Rr2-R.Rd2.*R.RL;
  R.Rr1 = sqrt(R.Rh1.^2 + R.Rw1.^2);
  R.Rh2 = R.Rh1 .* R.Rr2 ./ R.Rr1;
  R.Rw2 = sqrt(R.Rr1.^2 - R.Rh2.^2);
  R.Rd1 = R.Rr1 ./ R.RR1;
  R.Rs2 = R.Rw1 ./ R.RL;
  R.Rs1 = R.Rw2 ./ R.RL;
  R.s1 = R.Rs2;
  R.L = [];
  RLi = [];
  W = [
    R.s1.^2 + R.d1.^2
    -(2*R.r1.*R.d1 + (R.s1.^2 + R.d1.^2).*R.R2)
    (R.r1.^2 + R.R2.*R.r1.*R.d1)*ones(1,length(R.RL))
    ];
  for i=1:length(R.RL)
    rts = roots(W(:,i));
    rts = rts(imag(rts)==0);
    rts = rts(rts > 0);
    R.L = [R.L rts'];
    RLi = [RLi i*ones(1,length(rts))];
  end
  R.RL = R.RL(RLi);
  R.Rh1 = R.Rh1(RLi);
  R.Rr1 = R.Rr1(RLi);
  R.Rh2 = R.Rh2(RLi);
  R.Rw2 = R.Rw2(RLi);
  R.Rd1 = R.Rd1(RLi);
  R.Rs2 = R.Rs2(RLi);
  R.Rs1 = R.Rs1(RLi);
  R.s1 = R.Rs2;

  if ~isempty(R.L)
    R.w2 = R.s1.*R.L;
    R.h2 = R.r1 - R.d1.*R.L;
    R.r2 = sqrt(R.w2.^2 + R.h2.^2);
    R.d2 = R.r2 ./ R.R2;
    R.h1 = R.r2 - R.d2.*R.L;
    if any(R.h1 >= R.r1 | R.h1 <= 0)
      vh = R.h1 < R.r1 & R.h1 > 0;
      R.L = R.L(vh);
      R.w2 = R.w2(vh);
      R.h2 = R.h2(vh);
      R.r2 = R.r2(vh);
      R.d2 = R.d2(vh);
      R.h1 = R.h1(vh);
      R.RL = R.RL(vh);
      R.Rh1 = R.Rh1(vh);
      R.Rr1 = R.Rr1(vh);
      R.Rh2 = R.Rh2(vh);
      R.Rw2 = R.Rw2(vh);
      R.Rd1 = R.Rd1(vh);
      R.Rs2 = R.Rs2(vh);
      R.Rs1 = R.Rs1(vh);
      R.s1 = R.Rs2;
    end
    R.w1 = sqrt(R.r1.^2 - R.h1.^2);
    R.s2 = R.w1 ./ R.R1;
    R.r15 = (R.s1 .* R.r1)/tan(deg2rad(15));
  end
elseif isfield(P,'r2')
  % This and exexparam2a which uses it are illconceived.
  R.h2 = R.r1.*(1 - R.L./R.R1);
  if length(R.r2) > 1
    R.r2 = R.r2(R.h2 < R.r2);
  elseif length(R.h2) > 1
    R.h2 = R.h2(R.h2 < R.r2);
  elseif R.h2 >= R.r2
    R.h2 = [];
  end
  if isempty(R.r2) || isempty(R.h2)
    return;
  end
  R.h1 = R.h2.*R.r1./R.r2;
  R.w1 = sqrt(R.r1.^2 - R.h1.^2);
  R.w2 = sqrt(R.r2.^2 - R.h2.^2);
  R.s1 = R.w2./R.L;
  R.r15 = (R.s1 .* R.r1)/tan(deg2rad(15));
  % R.C2 = (R.r2-R.h1)./(R.r2.*R.L);
  R.d1 = R.r1./R.R1;
  
  R.Rr2 = R.r1;
  R.Rd2 = -n*R.d1;
  R.Rs2 = R.s1;
  R.RR2 = -R.R1/n;
  R.RL = R.Rw1./R.Rs2;
  R.Rh1 = R.Rr2 - R.Rd2.*R.RL;
  R.Rr1 = sqrt(R.Rw1.^2 + R.Rh1.^2);
  R.Rh2 = R.Rh1.*R.Rr2./R.Rr1;
  R.RR1 = R.Rr1.*R.RL./(R.Rr1-R.Rh2);
  R.d2 = R.r2./R.R2;
  R.s2 = R.w1./R.L;
  R.Rd1 = R.Rr1./R.RR1;
  R.Rw2 = R.Rr2*sqrt(1-((R.RR1-R.RL)*(R.RR2-R.RL))/(R.RR1*R.RR2));
  R.Rs1 = R.Rw2./R.RL;
end
