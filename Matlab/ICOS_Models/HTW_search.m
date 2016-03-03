%%
% HTW_search.m
% Need to accomodate possible change in n for 3765-3780 cm-1
%   Somewhere in the range of 2.4465 to 2.4467 based on interpolation
rd = 0.08;
th = 13;
mnc = 'HTW'; % sprintf('HTW_rd%02d_th%.1f', floor(rd*100), th);
SR = ICOS_sr_search('mnc',mnc,'B',0.4,'Rw1',0.2, 'rd', rd, 'th', th, ...
  'R1',[75 100 150],'R2',[75 100 150],'L',[50 100],'n',2.4466);
SR.enumerate;
%%
SR.design_tolerance;
%%
SR.build_tolerance;
%%
SR.explore_build_tolerance; % BT vs RLmin
SR.explore_build_tolerance(2); % BT vs DT
SR.explore_build_tolerance(3); % BT vs r2/r1
% SR.explore_build_tolerance(4); % R1 vs R2 (uninteresting)
%%
SR.explore_build_tolerance(5); % L
%%
SR.savefile;
%%
SR.focus('focus_visible',0);
%%
results = dir('IS_HTW*.mat');
for i=1:length(results)
  load(results(i).name);
  fprintf(1,'%d: %d focus solutions: %s \n', i, length(IS.res2), results(i).name);
  fprintf(1,'%d: d2 = %.4f\n', i, IS.res1(1).d2);
  check_params(i, IS.res1(1));
  for j=1:length(IS.res2)
    fprintf(1, '  %d: ', j);
    Lenses = IS.res2(j).Lenses;
    for k=1:length(Lenses)
      fprintf(1, ' %s', Lenses{k});
    end
    fprintf(1,'\n');
  end
end
%%
results = dir('IS_HTW*.mat');
for i=1:length(results)
  load(results(i).name);
  fprintf(1,'%d: RL: %.1f R1: %.0f R2: %.0f L: %.1f r2/r1: %.1f\n', i, ...
    IS.res1.RL, IS.res1.R1, IS.res1.R2, IS.res1.L, IS.res1.r2/IS.res1.r1);
end
