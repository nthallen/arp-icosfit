%%
% scratch_151216.m
%  Look at different focuses
SolScale = 4;
oldfile = sprintf('IS_L50b_31.2.%d.1_31.2.mat', SolScale);
load(oldfile);
newmnc = sprintf('L50c.%d', SolScale);
IS.ISopt.mnc = newmnc;
%%
IS.explore_focus;
%%
IS.savefile
%%
rng_state = rng;
% If analyze has to be interrupted, we'll need to retrieve rng_state
% from one of the early results
IS.analyze('rng_state', rng_state);

%%
res = collect_results('files','IS_L50c*.mat');
%%
Lrat = zeros(length(res),1);
idx = [res.index];
for i=1:length(res)
  LS = IS.res2(idx(i)).Lens_Space(end);
  DS = IS.res2(idx(i)).detector_spacing;
  Lrat(i) = LS/(LS+DS);
end
max_pwr = [res.max_pwr];
plot(Lrat, max_pwr, '*');
  
