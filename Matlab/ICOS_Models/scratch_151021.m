%% scratch_151021.m
% Look at detector/preamp bandwidth as alternate constraint from coherence
% length.
%
% Etalon produces fringes with spacing of c/(2nL). We'll assume n=1.
TR = 1000; % tuning rate in cm-1/sec
BW = 1e6; % detector/preamp bandwidth, Hz
c = lightspeed; % cm/sec
fsrmax = TR/BW; % cm-1
Lmin = 1/(2*fsrmax);
