%%
% The basic idea here is:
%   Given R1, R2, RR1 and two other parameters, we should be
%   able to solve for the entire RIM/ICOS configuration.
%   Since R1, R2 and RR1 come in discrete values, it makes sense
%   to fix those. For the two additional values to fix, I think
%   it makes sense to fix r1 and Rw1. Generally configurations will
%   scale with r1, so this just establishes the scale of the solution.
%   Rw1 defines the minimum beam spacing to make the Herriott mirror
%   work. Since RL = Rw1/Rs2, it's clear that for a given ICOS
%   configuration, a larger Rw1 will give us a longer RIM, and we'd
%   like to keep it at a minimum.
%% Radius of curvature for plano-convex 1" optics
% These are some of the possible values for R2
ZCPXR = -0.1 * [
  71.46
  89.14
  106.15
  139.96
  180.32
  209.91
  281.20
  356.51
  698.17
];
% ZCPXC = [1./ZCPXR; 0];
IR2 = ZCPXR;
HFL1 = [ % in inches
  4
  6
  8
  12
  18
  20
  24
  30
  36
  40
  ];
HR1 = HFL1 * 2 / 2.54;
nIR2 = length(IR2);
nHR1 = length(HR1);
ntrials = nIR2*nHR1;

P.R1 = 75;
P.r1 = 3;
P.Rw1 = 0.2;
res(ntrials) = exparam(P); % should be empty, just to initialize
trialn = 1;
for IR2i = 1:nIR2
  P.R2 = IR2(IR2i);
  for HR1i = 1:nHR1
    P.RR1 = HR1(HR1i);
    res(trialn) = exparam(P);
    trialn = trialn+1;
  end
end

vL = arrayfun(@(x) length(x{1}), {res.L}) > 0;
res = res(vL);
