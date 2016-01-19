%%
% scratch_151222.m ICOS diagram figure for paper
% Try with L = 20, short coherence length
SR = ICOS_sr_search('mnc','L20','L',20,'B',0.4,'Rw1',0.2,'C',1000);
SR.enumerate;
%%
SR.design_tolerance;
SR.build_tolerance;
%%
SR.explore_build_tolerance(3);
%%
SR.savefile;
%%
SR.focus('focus_visible',0);
%%
load('IS_L20.103_36.1.mat');
IS.explore_focus;
%%
n = find([IS.res2.sel]);
%%
P = render_model(IS.res2(n(1)),'visibility',0,'visible',1,'HR',0, ...
  'view',[0 90],'ICOS_passes_per_injection', 30);
PM = ICOS_Model6(P);
%%
% Now format it for publication
xl = xlim;
xl = [-3 32];
xlim(xl);
ylim(P.r1 * [-1.1 1.1]);
%% OK, let's try that again. I want a cylindrical diagram for this
set(gca,'visible','off');
set(gcf,'color',[1 1 1]);
%%
xl = xlim;
yl = ylim;
set(gca,'units','pixels');
set(gcf,'units','pixels');
AP = [0 0 1000 0];
FP = get(gcf,'position');
AP2 = [0, 0, AP(3), AP(3)*diff(yl)/diff(xl)];
FP2 = [FP(1:2) AP2(3:4)];
set(gca,'position',AP2);
set(gcf,'position',FP2);
shg;
%%
PP = get(gcf,'PaperPosition');
PP(4) = PP(3)*diff(yl)/diff(xl);
set(gcf,'PaperPosition',PP);
%%
print('figure_cyl_ICOS.png','-dpng','-r300');

%%
P = render_model(IS.res2(n(1)),'visibility',[],'visible',1, ...
  'view',[0 90],'ICOS_passes_per_injection', 1);
PM = ICOS_Model6(P);
%%
% Now format it for publication
xl = xlim;
xl = [-23 32];
xlim(xl);
ylim(P.r1 * [-1.1 1.1]);
%% OK, let's try that again. I want a cylindrical diagram for this
set(gca,'visible','off');
set(gcf,'color',[1 1 1]);
%%
xl = xlim;
yl = ylim;
set(gca,'units','pixels');
set(gcf,'units','pixels');
AP = [0 0 1000 0];
FP = get(gcf,'position');
AP2 = [0, 0, AP(3), AP(3)*diff(yl)/diff(xl)];
FP2 = [FP(1:2) AP2(3:4)];
set(gca,'position',AP2);
set(gcf,'position',FP2);
shg;
%%
PP = get(gcf,'PaperPosition');
PP(4) = PP(3)*diff(yl)/diff(xl);
set(gcf,'PaperPosition',PP);
%%
print('figure_cyl_RIM.png','-dpng','-r300');

