%% scratch_160127: Model of benchtop HCl/ME system with crappy mirrors
% Create figure similar to Figure 17 (overfocusing)
% We can't use ICOS_sr_search because it will just fail, but we reimplement
% some of the core processing here.
R = 75;
L = 50;
w_r = sqrt(L*(2*R-L)/R^2);
phi = asind(w_r);
r1 = 0.75*2.54; % 1" through Lens1 is overfocused
w1 = r1 * w_r;
B = 0.4;
Rw1 = 1;
s1 = w1/L;
n = 2.4361;
Rs2 = s1;
RL = Rw1/Rs2;
SP = struct('R1',R,'R2',R,'L',L,'RL',RL,'Rw1',Rw1);
Res = exparam(SP);
%%
% IS = ICOS_search('R1',R,'R2',R,'L',L,
IS = ICOS_search('mnc', 'JWa','R1',SP.R1,'R2',SP.R2,'L',SP.L, ...
  'RR1',Res.RR1,'Rw1',SP.Rw1, 'RL_lim', [0.90,1.1]*SP.RL);
%%
IS.search_ICOS_RIM;
%%
load IS_JWa.mat
%%
% This did not return any successful focuses
% fix = {'Lens1','ZC_PM_12_12' };
% for th=[12 13 14 15 16]
%   IS.search_focus2('max_lenses',2,'fix_lenses',fix,'det_acc_limit',th);
% end
fix = {'Lens1','ZC_PM_12_25' };
% ths = [14.8 15 16];
ths = 13:17;
for th=ths
  IS.search_focus2('max_lenses',2,'fix_lenses',fix,'det_acc_limit',th,...
    'allow_nondecreasing_focus',1);
end
%%
% Let's look at what's going on:
close all
P = render_model(IS.res1);
P.Lenses = {'Lens1', 'ZC_PM_12_25'};
P.Lens_Space = [0.1 8.6];
P.HR = 0;
P.focus = 1;
P.detector_spacing = 1;
P.herriott_spacing = 3;
P.evaluate_endpoints = 0;
P.visibility = [0 0 0];
P.view = [0 90];
figure;
PM = ICOS_Model6(P);
%
[xyz, dxyz, d, s] = PM.M.extract_endpoints_skew(4 + length(P.Lenses));
d_angle = atand(sqrt(d.^2 + s.^2));
figure;
plot(d_angle,'.');
shg;
%%
% We are already overfocused, coming in at 16.4 deg
fix = {'Lens1'};
for th=[16.5 17]
  IS.search_focus2('max_lenses',2,'fix_lenses',fix,'det_acc_limit',th, ...
    'allow_nondecreasing_focus',1);
end

%%
% Create a figure
% First do the analysis
IS.explore_focus;
IS.savefile;
IS.analyze('HR',0,'T',.005);
%%
close all;
clear;
load IS_JWa.mat
%%
% Show effects of focusing including lens diagram. These samples differ only on focal angle
V = 1:5;
Naxes = length(V);
Opt.ICOS_passes = 50;
Opt.Nsamples = 100;
CL = [0 0];
XL = [];
YL = [];
tmpfig = cell(Naxes,2);
for ni = 1:Naxes
  tmpfig{ni,1} = sprintf('tmp_JWa_%d.ray.fig',ni);
  tmpfig{ni,2} = sprintf('tmp_JWa_%d.heat.fig',ni);
  i = V(ni);
  IBmnc = sprintf('%s.%d_%dx%d', IS.ISopt.mnc, ...
    IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
  ofile = sprintf('IB_%s.mat', IBmnc);
  load(ofile);
  ff = IB.Integrate;
  delete(ff([1 3 4]));
  fg = ff(2);
  ax = findobj(fg,'type','axes');
  colorbar(ax,'off');
  title(ax,'');
  xlabel(ax,'');
  ylabel(ax,'');
  set(ax,'units','pixels');
  cl = get(ax,'Clim');
  CL(2) = max(cl(2), CL(2));
  drawnow;
  savefig(fg,tmpfig{ni,2});
  delete(fg);
  %
  % ICOS injections selected here by trial and error
  fg = figure;
  IB.draw('ICOS_passes_per_injection',23,'visibility',[0 0 0],'HR',0,'visible',1,'view',[0 90]);
  ax = findobj(fg,'type','axes');
  xl = get(ax,'xlim');
  yl = get(ax,'ylim');
  if isempty(XL)
    XL = xl;
  else
    XL(1) = min(XL(1),xl(1));
    XL(2) = max(XL(2),xl(2));
  end
  if isempty(YL)
    YL = yl;
  else
    YL(1) = min(YL(1),yl(1));
    YL(2) = max(YL(2),yl(2));
  end
  drawnow;
  savefig(fg,tmpfig{ni,1});
  delete(fg);
  %
  delete(IB);
end
save tmp_JWa.mat XL YL CL tmpfig V IS
clear
%%
clear
close all
load tmp_JWa.mat
%
% #### Figure out how to arrange these
% Will set the heatmaps to the same size as the ray trace
% Need to figure out the true size of the ray trace plot box.
% They should now all be the same size
newfig = figure;
menuheight = 80;
taskbarheight = 40;
figborder = 9;
topspace = 20;
bottomspace = 70;
leftspace = 20;
midspace = 10;
rightspace = 20;
Naxes = length(V);
screen = get(groot,'screensize');
axesheight = screen(4)-menuheight-taskbarheight - topspace - bottomspace;
rayheight = axesheight/Naxes;
raywidth = rayheight * diff(XL)/diff(YL);
figwidth = leftspace+raywidth+midspace+rayheight+rightspace;
figheight = topspace+Naxes*rayheight+bottomspace;
figpos = [figborder,screen(4)-figheight-menuheight,figwidth,figheight];
set(newfig,'position',figpos,'color',[1 1 1]);
drawnow;
get(newfig,'position')
%
Y = 0.1*[-1 1 1 -1 -1 ];
Z = 0.1*[-1 -1 1 1 -1 ];
y = figheight-topspace;
naxes = zeros(Naxes,1);
nfxes = zeros(Naxes,1);
ni = 1;
%
for ni=1:Naxes
  %%
  y = y-rayheight;
  fig = openfig(tmpfig{ni,1});
  ax = findobj(fig,'type','axes');
  nfx = copyobj(ax,newfig); % this is the ray trace
  delete(fig);
  set(nfx,'Xlim',XL,'Ylim',YL,'units','pixels');
  set(nfx,'position',[leftspace,y,raywidth,rayheight]);
  set(nfx,'visible','off');
  %
  fig = openfig(tmpfig{ni,2});
  ax = findobj(fig,'type','axes');
  nax = copyobj(ax,newfig); % this is the heat map
  delete(fig);
  %
  set(nax,'Clim',CL);
  naxes(ni) = nax;
  nfxes(ni) = nfx;
  %
  set(nax,'units','pixels');
  set(nax,'position',[leftspace+raywidth+midspace,y,rayheight,rayheight], ...
    'visible','on','XTickLabel',[],'YTickLabel',[]);
  hold(nax,'on');
  plot(nax,Y,Z,'w');
  hold(nax,'off');
  h = ylabel(nax,sprintf('%.0f^{\\circ}',IS.res2(V(ni)).theta));
  set(h,'horizontalalignment','right','verticalalignment','middle', ...
    'rotation',0,'fontsize',14,'fontweight','bold');
end
drawnow;
%
pos = get(nax,'position');
H = colorbar('peer',nax,'southoutside');
drawnow;
set(nax,'position',pos);
drawnow;
% set(H,'units','pixels');
% pos = get(H,'position');
% pos(2) = pos(2)+25;
% set(H,'position',pos);
set(H,'fontsize',10,'fontweight','bold');
%%
print_figure_v1(newfig,'figure_JWa.png',[8 10.5],[naxes;H]);
%%
clear
%%
load IS_JWa.mat
P = render_model(IS.res1,'HR',0,'herriott_spacing',8);
%%
draw_targets(P,1/2.54); % Scale to match page units in inches
%%
print -dpng -r300 target_JWa.png

%% How to examine optic positioning
load IB_JWa.1_50x100.mat
PM = IB.draw('ICOS_passes_per_injection',23,'visibility',0, ...
  'HR',0,'visible',1,'view',[0 90],'herriott_spacing',8);
%% How to examine optic positioning
load IB_JWa.5_50x100.mat
PM = IB.draw('ICOS_passes_per_injection',23,'visibility',0, ...
  'HR',0,'visible',1,'view',[0 90],'herriott_spacing',8);
