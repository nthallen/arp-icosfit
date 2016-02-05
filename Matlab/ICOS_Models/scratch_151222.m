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
%%
IS.explore_focus;
%%
IS.savefile;
%%
n = find([IS.res2.sel]);
%%
P = render_model(IS.res2(n(1)),'visibility',0,'visible',1,'HR',0, ...
  'view',[0 90],'ICOS_passes_per_injection', 30, 'herriott_spacing', 3);
PM = ICOS_Model6(P);
%% OK, let's try that again. I want a cylindrical diagram for this
set(gca,'visible','off');
set(gcf,'color',[1 1 1]);
% Now do annotations
% Draw laser box
sth = P.dy/sqrt(P.dy^2+1);
cth = 1/sqrt(P.dy^2+1);
wtdl = 0.5;
ltdl = 1.5;
ytdl = (wtdl/2)*[-1 1 1 -1 -1];
xtdl = ltdl * [0 0 -1 -1 0];
rtdl = [xtdl' ytdl']*[cth -sth; sth cth];
ttdl = [rtdl(:,1)-P.herriott_spacing-P.CT1, rtdl(:,2)+P.y0+P.dy*P.herriott_spacing];
hold on;
fill(ttdl(:,1),ttdl(:,2),[.5 .5 1]);
plot(ttdl(:,1),ttdl(:,2),'k');
hold off;
dt = .4;
t1 = text(-P.herriott_spacing-P.CT1-ltdl/2, P.y0-dt, 'TDL');
set(t1,'horizontalalignment','center','verticalalignment','top', ...
  'fontsize',14,'fontweight','bold');
%
t2 = text(dt, P.r1, 'M_1');
set(t2,'horizontalalignment','left','verticalalignment','top', ...
  'fontsize',14,'fontweight','bold');
%
t3 = text(P.mirror_spacing-dt, P.r2, 'M_2');
set(t3,'horizontalalignment','right','verticalalignment','middle', ...
  'fontsize',14,'fontweight','bold');
%
Lens = P.LensTypes.(P.Lenses{1});
tx = P.mirror_spacing + P.CT1 + P.Lens_Space(1)+Lens.CT;
Lens = P.LensTypes.(P.Lenses{1});
t4 = text(tx, Lens.r+dt/2, 'L_1');
set(t4,'horizontalalignment','center','verticalalignment','bottom', ...
  'fontsize',14,'fontweight','bold');
%
Lens = P.LensTypes.(P.Lenses{2});
tx = tx + P.Lens_Space(2) + Lens.CT;
t5 = text(tx, Lens.r+dt/2, 'L_2');
set(t5,'horizontalalignment','center','verticalalignment','bottom', ...
  'fontsize',14,'fontweight','bold');
%
tx = tx + P.detector_spacing;
t6 = text(tx, P.D_l/2 + dt, 'D');
set(t6,'horizontalalignment','Center','verticalalignment','bottom', ...
  'fontsize',14,'fontweight','bold');
%%
% Now format it for publication
xl = xlim;
xl = [min(ttdl(:,1))-2*dt tx+2*dt];
%xl = [-3 32];
xlim(xl);
ylim(P.r1 * [-1.2 1.2]);
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
print_figure_v1(gca, 'figure_cyl_ICOS.png', [8 5], [t1 t2 t3 t4 t5 t6]);
% %%
% PP = get(gcf,'PaperPosition');
% PP(4) = PP(3)*diff(yl)/diff(xl);
% set(gcf,'PaperPosition',PP);
% %%
% print('figure_cyl_ICOS.png','-dpng','-r300');

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

