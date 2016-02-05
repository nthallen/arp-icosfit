%% scratch_160120.m
% Need to reanalyze to see how acceptance angle affects choices
rd = 0.08;
th = 13;
mnc = sprintf('L50d_rd%02d_th%.1f', floor(rd*100), th);
SR = ICOS_sr_search('mnc',mnc,'B',0.4,'Rw1',0.2, 'rd', rd, 'th', th);
SR.enumerate;
SR.design_tolerance;
SR.build_tolerance;
SR.explore_build_tolerance;
SR.explore_build_tolerance(2);
SR.explore_build_tolerance(3);
SR.explore_build_tolerance(4);
%%
SR.savefile;
%%
SR.focus('focus_visible',0);
%%
files = dir('IS_L50d_rd08_*.mat');
for i=1:length(files)
  clear IS
  load(files(i).name);
  if exist('IS','var') && ~isempty(IS.res2)
    IS.explore_focus;
    IS.savefile;
  end
end
%%
files = dir('IS_L50d_rd08*.mat');
for i=1:length(files)
  clear IS
  load(files(i).name);
  if exist('IS','var') && ~isempty(IS.res2) && any([IS.res2.sel])
    fprintf(1,'Analyzing %s\n', IS.ISopt.mnc);
    IS.analyze;
    close all
  end
end
%%
res = collect_results('files', 'IS_L50d_rd08*.mat');
%%
% Select the best option
resn = 2;
mnc = res(resn).mnc;
%%
% Create graphics for full integration
pat = sprintf('IB_%s.%d*.mat', mnc, res(resn).index);
files = dir(pat);
load(files(1).name);
ff = IB.Integrate;
%%
delete(ff(3:4));
ff = ff(1:2);
%%
Y = 0.1*[-1 1 1 -1 -1];
Z = 0.1*[-1 -1 1 1 -1];
ax = findobj(ff(1),'type','axes');
hold(ax,'on');
plot(ax,Y,Z,'w');
hold(ax,'off');
set(ax,'fontsize',12,'fontweight','bold');
title(ax,'High Resolution');
xlabel(ax,'Y cm');
ylabel(ax,'Z cm');
print_figure_v1(ax, 'figure_HR_focused.png', [3 3], ax);
%%
ax = findobj(ff(2),'type','axes');
set(ax,'fontsize',12,'fontweight','bold');
title(ax,'Detector Resolution');
xlabel(ax,'Y cm');
ylabel(ax,'Z cm');
print_figure_v1(ax, 'figure_DR_focused.png', [3 3], ax);
%%
P = IB.P;
P.HR = 0;
P.visibility = [0];
P.herriott_spacing = 3;
%%
% PM = ICOS_Model6(P);
%%
Opt.Nsamples = 100;
Opt.ICOS_passes = 50;
IBmnc = sprintf('%s.NR_%dx%d', mnc, ...
  Opt.ICOS_passes, Opt.Nsamples);
IB = ICOS_beam(@ICOS_Model6,P);
%%
IB.Sample('beam_samples', Opt.Nsamples, ...
  'ICOS_passes', Opt.ICOS_passes, 'opt_n', 6, ...
  'n_optics', 6, 'Track_Power', 1, ...
  'mnc', IBmnc);
ff = IB.Integrate;
IB.savefile;
%%
delete(ff(3:4));
ff = ff(1:2);
%%
Y = 0.1*[-1 1 1 -1 -1];
Z = 0.1*[-1 -1 1 1 -1];
ax = findobj(ff(1),'type','axes');
hold(ax,'on');
plot(ax,Y,Z,'w');
hold(ax,'off');
set(ax,'fontsize',12,'fontweight','bold');
title(ax,'High Resolution');
xlabel(ax,'Y cm');
ylabel(ax,'Z cm');
print_figure_v1(ax, 'figure_HR_NR_focused.png', [3 3], ax);
%%
ax = findobj(ff(2),'type','axes');
set(ax,'fontsize',12,'fontweight','bold');
title(ax,'Detector Resolution');
xlabel(ax,'Y cm');
ylabel(ax,'Z cm');
print_figure_v1(ax, 'figure_DR_NR_focused.png', [3 3], ax);
