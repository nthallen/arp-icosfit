function R = exparam(P)
% R = exparam(P)
% P is a struct defining R1, r1, r2 and L (and Rw1).
% or P defines R1, R2, RR1, r1 and Rw1
% One of these may be a vector.
Res = cell(1,0);
R.n = 2.4361;

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

if isfield(P,'RR1') && isfield(P,'r1') % assume RR1, R1, R2, r1 and Rw1
  R.d1 = R.r1./R.R1;
  R.Rd2 = -R.d1*R.n;
  R.RR2 = - R.R1/R.n;
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
  for RLi = 1:length(rts)
    R.RL = rts(RLi);
    
    R.Rh1 = R.Rr2-R.Rd2.*R.RL;
    R.Rr1 = sqrt(R.Rh1.^2 + R.Rw1.^2);
    R.Rh2 = R.Rh1 .* R.Rr2 ./ R.Rr1;
    R.Rw2 = sqrt(R.Rr2.^2 - R.Rh2.^2);
    R.Rd1 = R.Rr1 ./ R.RR1;
    R.Rs2 = R.Rw1 ./ R.RL;
    R.Rs1 = R.Rw2 ./ R.RL;
    R.s1 = R.Rs2;
    W = [
      R.s1.^2 + R.d1.^2
      -(2*R.r1.*R.d1 + (R.s1.^2 + R.d1.^2).*R.R2)
      (R.r1.^2 + R.R2.*R.r1.*R.d1)
      ];
    Lrts = roots(W);
    Lrts = Lrts(imag(Lrts)==0);
    Lrts = Lrts(Lrts > 0);
    for Li = 1:length(Lrts)
      R.L = Lrts(Li);
      R.s1 = R.Rs2;
      R.w2 = R.s1.*R.L;
      R.h2 = R.r1 - R.d1.*R.L;
      R.r2 = sqrt(R.w2.^2 + R.h2.^2);
      R.d2 = R.r2 ./ R.R2;
      R.h1 = R.r2 - R.d2.*R.L;
      if R.h1 < R.r1 && R.h1 > 0
        R.w1 = sqrt(R.r1.^2 - R.h1.^2);
        R.s2 = R.w1 ./ R.L;
        R.r15 = (R.s1 .* R.r1)/tan(deg2rad(15));
        Res{1,end+1} = R;
      end
    end
  end
elseif isfield(P,'RR1') && isfield(P,'L')
  % assuming RR1, R1, R2, Rw1, L
  if (R.R1+R.R2-R.L)/(R.R1*R.R2) > 0
    R.RR2 = -R.R1/R.n;
    k = sqrt((R.R1-R.L)*(R.R1+R.R2-R.L)/(R.L*R.R1^2*(R.R2-R.L)));
    V = [
      R.Rw1*(R.RR1+R.RR2)
      -R.Rw1^2*((2/k) + (R.RR1/(k*R.RR2)) + k*R.RR1*R.RR2)
      R.Rw1^3*((1/(k^2*R.RR2))+R.RR2)
      ];
    rts = roots(V);
    rts = rts(imag(rts)==0);
    rts = rts(rts > 0);
    for r1i = 1:length(rts)
      R.r1 = rts(r1i);
      R.r2 = R.r1*sqrt((R.R2*(R.R1-R.L))/(R.R1*(R.R2-R.L)));
      R.h1 = R.r1*sqrt((R.R1-R.L)*(R.R2-R.L)/(R.R1*R.R2));
      if R.h1 < R.r1
        R.w1 = sqrt(R.r1^2-R.h1^2);
        R.h2 = R.h1*R.r2/R.r1;
        R.w2 = R.w1*R.r2/R.r1;
        R.d1 = R.r1/R.R1;
        R.d2 = R.r2/R.R2;
        R.s1 = R.w2/R.L;
        R.s2 = R.w1/R.L;
        R.Rd2 = -R.n*R.d1;
        R.Rs2 = R.s1;
        R.Rr2 = R.r1;
        R.RL = R.Rw1/R.Rs2;
        if (R.RR1+R.RR2-R.RL)/(R.RR1*R.RR2) > 0
          %R.Rr1 = R.Rr2*sqrt((R.RR1*(R.RR2-R.RL))/(R.RR2*(R.RR1-R.RL)));
          R.Rh1 = R.Rr2*(1-R.RL/R.RR2);
          R.Rr1 = sqrt(R.Rw1^2 + R.Rh1^2);
          %R.Rh2 = R.Rh1*R.Rr2/R.Rr1;
          R.Rh2 = R.Rr1*(1 - R.RL/R.RR1);
          R.Rw2 = R.Rw1*R.Rr2/R.Rr1;
          R.Rd1 = R.Rr1/R.RR1;
          R.Rs1 = R.Rw2/R.RL;
          R.r15 = (R.s1 .* R.r1)/tan(deg2rad(15));
          Res{1,end+1} = R;
        end
      end
    end
  end
end
R = [Res{:}]';
