%% scratch_160113: Model of benchtop HCl/ME system with crappy mirrors
% Includes figure for focusing
% We can't use ICOS_sr_search because it will just fail, but we reimplement
% some of the core processing here.
R = 75;
L = 50;
w_r = sqrt(L*(2*R-L)/R^2);
phi = asind(w_r);
r1 = 1*2.54;
w1 = r1 * w_r;
B = 0.4;
Rw1 = 1;
s1 = w1/L;
n = 2.4361;
Rs2 = s1;
RL = Rw1/Rs2;
SP = struct('R1',R,'R2',R,'L',L,'RL',RL,'Rw1',Rw1);
Res = exparam(SP);
% IS = ICOS_search('R1',R,'R2',R,'L',L,
IS = ICOS_search('mnc', 'JW','R1',SP.R1,'R2',SP.R2,'L',SP.L, ...
  'RR1',Res.RR1,'Rw1',SP.Rw1, 'RL_lim', [0.90,1.1]*SP.RL);
%%
IS.search_ICOS_RIM;
%%
IS.search_focus2('max_lenses',2);
%% Do analysis with bad mirrors, 50 passes
V = find([IS.res2.sel]);
ff = [];
for ni = 1:length(V)
  i = V(ni);
  P = render_model(IS.res2(i));
  P.HR = 0; % No RIM for now
  P.herriott_spacing = 3;
  P.T = .005; % 5000 ppm
  Opt.ICOS_passes = 50;
  Opt.Nsamples = 100;
  IBopt = {};
  
  n_optics = 4 + length(P.Lenses);
  opt_n = n_optics;
  Track_Power = 1;
  IBmnc = sprintf('%s.LowR.%d_%dx%d', IS.ISopt.mnc, ...
    IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
  ofile = sprintf('IB_%s.mat', IBmnc);
  if ~exist(ofile, 'file')
    P.visible = 0;
    % P.HR = 0;
    P.focus = 1;
    P.evaluate_endpoints = -1;
    IB = ICOS_beam(@ICOS_Model6, P);
    IB.Sample('beam_samples', Opt.Nsamples, ...
      'ICOS_passes', Opt.ICOS_passes, 'opt_n', opt_n, ...
      'n_optics', n_optics, 'Track_Power', Track_Power, ...
      'mnc', IBmnc, IBopt{:});
    if opt_n == n_optics
      ff2 = IB.Integrate;
      if ~isempty(ff)
        delete(ff);
        drawnow;
      end
      ff = ff2;
    end
    IB.savefile;
    % save(ofile, 'IB');
    % fprintf(1, 'ICOS_beam %d Saved result to %s\n', i, ofile);
    Pwr = IB.PowerSummary;
    IS.res2(i).NH = Pwr.NH;
    if isfield(Pwr,'max_pwr')
      IS.res2(i).max_pwr = Pwr.max_pwr;
    elseif ~isfield(IS.res2,'max_pwr')
      IS.res2(i).max_pwr = [];
    end
    IS.savefile;
  else
    fprintf(1, 'Skipping result %d\n', i);
  end
end
%% Create some additional focuses
fix = IS.res2(17).Lenses;
for th=[14 16 17]
  IS.search_focus2('max_lenses',2,'fix_lenses',fix,'det_acc_limit',th);
end
%% Do analysis with good mirrors, 50 passes, no RIM
% V = find([IS.res2.sel]);
V = 19:21;
V = 17;
V = [18 19 17 20 21];
ff = [];
for ni = 1:length(V)
  i = V(ni);
  P = render_model(IS.res2(i));
  P.HR = 0; % No RIM for now
  P.herriott_spacing = 3;
  % P.T = .005; % 5000 ppm
  Opt.ICOS_passes = 1000;
  Opt.Nsamples = 100;
  IBopt = {'Int', struct('Yr',[-.72 .72;-.72 .72],'Zr',[-.72 .72; -.72 .72])};
  
  n_optics = 4 + length(P.Lenses);
  opt_n = n_optics;
  Track_Power = 1;
  IBmnc = sprintf('%s.HiR.%d_%dx%d', IS.ISopt.mnc, ...
    IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
  ofile = sprintf('IB_%s.mat', IBmnc);
  if ~exist(ofile, 'file')
    P.visible = 0;
    % P.HR = 0;
    P.focus = 1;
    P.evaluate_endpoints = -1;
    IB = ICOS_beam(@ICOS_Model6, P);
    IB.Sample('beam_samples', Opt.Nsamples, ...
      'ICOS_passes', Opt.ICOS_passes, 'opt_n', opt_n, ...
      'n_optics', n_optics, 'Track_Power', Track_Power, ...
      'mnc', IBmnc, IBopt{:});
    if opt_n == n_optics
      ff2 = IB.Integrate;
      if ~isempty(ff)
        delete(ff);
        drawnow;
      end
      ff = ff2;
    end
    IB.savefile;
    % save(ofile, 'IB');
    % fprintf(1, 'ICOS_beam %d Saved result to %s\n', i, ofile);
    Pwr = IB.PowerSummary;
    IS.res2(i).NH = Pwr.NH;
    if isfield(Pwr,'max_pwr')
      IS.res2(i).max_pwr = Pwr.max_pwr;
    elseif ~isfield(IS.res2,'max_pwr')
      IS.res2(i).max_pwr = [];
    end
    IS.savefile;
  else
    fprintf(1, 'Skipping result %d\n', i);
  end
end
%%
close all
%%
load IS_JW.mat
%%
% Show effects of focusing. These samples differ only on focal angle
V = [18 19 17 20 21];
Opt.ICOS_passes = 1000;
Opt.Nsamples = 100;
fg = zeros(length(V),1);
ax = zeros(length(V),1);
CL = [0 0];
for ni = 1:length(V)
  i = V(ni);
  IBmnc = sprintf('%s.HiR.%d_%dx%d', IS.ISopt.mnc, ...
    IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
  ofile = sprintf('IB_%s.mat', IBmnc);
  load(ofile);
  ff = IB.Integrate;
  delete(ff([1 3 4]));
  fg(ni) = ff(2);
  ax(ni) = findobj(gcf,'type','axes');
  cl = get(ax(ni),'Clim');
  CL(2) = max(cl(2), CL(2));
end
set(ax,'Clim',CL);
%% Clean up the figures we've saved
for ni=1:length(ax)
  xl = get(ax(ni),'XLim');
  yl = get(ax(ni),'YLim');
  colorbar(ax(ni),'off');
  title(ax(ni),'');
  xlabel(ax(ni),'');
  ylabel(ax(ni),'');
  drawnow;
  set(ax(ni),'units','pixels');
%   pos = get(ax(ni),'position');
%   DAR = diff(yl)/diff(xl);
%   PAR = pos(4)/pos(3);
%   if PAR > DAR
%     dh = pos(4) - pos(3)*DAR;
%     fprintf(1,'%d: axes are %.1f pixels taller than expected\n', ni, dh);
%   else
%     dw = pos(3) - pos(4)/DAR;
%     fprintf(1,'%d: axes are %.1f pixels wider than expected\n', ni, dw);
%   end
end
%%
newfig = figure;
pos = get(newfig,'position');
pos(1) = 0;
pos(3) = 4600;
set(newfig,'position',pos);
drawnow;
pos = get(newfig,'position');
figwid = pos(3);
leftspace = 50;
rightspace = 120;
Naxes = length(ax);
axwid = (figwid-leftspace-rightspace)/Naxes;

x = leftspace;
Y = 0.1*[-1 1 1 -1 -1 ];
Z = 0.1*[-1 -1 1 1 -1 ];
naxes = zeros(length(ax),1);
for ni=1:length(ax)
  nax = copyobj(ax(ni),newfig);
  naxes(ni) = nax;
  pos = get(nax,'position');
  pos(1) = x;
  pos(3) = axwid;
  pos(4) = axwid;
  set(nax,'position',pos);
  if ni > 1
    set(nax,'YTickLabel',[]);
  end
  set(nax,'XTickLabel',[],'fontsize',14,'fontweight','bold');
  title(nax,sprintf('%.0f^{\\circ}',IS.res2(V(ni)).theta));
  hold(nax,'on');
  plot(nax,Y,Z,'w');
  hold(nax,'off');
  x = x + axwid;
end
drawnow;
pos = get(nax,'position');
colorbar('peer',nax);
drawnow;
set(nax,'position',pos);
%pos = get(newfig,'position');
%pos(3) = x+dw+30;
%set(newfig,'position',pos);
%%
print_figure_v1(newfig,'figure_overfocus.png',[8 4],naxes);
%%
close all
%%
clear;
load IS_JW.mat
%%
% Show effects of focusing including lens diagram. These samples differ only on focal angle
V = [18 19 17 20 21];
% V = 18;
Naxes = length(V);
Opt.ICOS_passes = 1000;
Opt.Nsamples = 100;
CL = [0 0];
XL = [];
YL = [];
tmpfig = cell(Naxes,2);
for ni = 1:Naxes
  tmpfig{ni,1} = sprintf('tmp_overfocus_%d.ray.fig',ni);
  tmpfig{ni,2} = sprintf('tmp_overfocus_%d.heat.fig',ni);
  i = V(ni);
  IBmnc = sprintf('%s.HiR.%d_%dx%d', IS.ISopt.mnc, ...
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
  IB.draw('ICOS_passes_per_injection',23,'visibility',[0 0 0],'visible',1,'view',[0 90]);
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
save tmp_overfocus.mat XL YL CL tmpfig V IS
clear
%%
clear
close all
load tmp_overfocus.mat
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
%%
set(H,'fontsize',10,'fontweight','bold');
%%
print_figure_v1(newfig,'figure_overfocus2.png',[8 10.5],[naxes;H]);

%%
% Produce Hi Res unfocused No RIM integration figures
clear
load IB_JW.HiR.17_1000x100.mat
%%
ff = IB.Integrate;
delete(ff(3:4));
ff = ff(1:2);
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
print_figure_v1(ax, 'figure_HR_NR_unfocused.png', [3 3], ax);
%%
ax = findobj(ff(2),'type','axes');
set(ax,'fontsize',12,'fontweight','bold');
title(ax,'Detector Resolution');
xlabel(ax,'Y cm');
ylabel(ax,'Z cm');
print_figure_v1(ax, 'figure_DR_NR_unfocused.png', [3 3], ax);


%% Focused example, no RIM
% This has been redone in scratch_160120.m as L50d
load('IS_L50c.4.mat');
i = 48; % by inspection...
ff = [];
  P = render_model(IS.res2(i));
  % P.HR = 0; % No RIM for now
  % P.T = .005; % 5000 ppm
  Opt.ICOS_passes = 100; % 100 for debugging
  Opt.Nsamples = 100;
  IBopt = {};
  
  n_optics = 4 + length(P.Lenses);
  opt_n = n_optics;
  Track_Power = 1;
  IBmnc = sprintf('%s.%d_%dx%d', IS.ISopt.mnc, ...
    IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
  ofile = sprintf('IB_%s.mat', IBmnc);
  if ~exist(ofile, 'file')
    P.visible = 0;
    P.HR = 0; % No RIM
    P.focus = 1;
    P.evaluate_endpoints = -1;
    IB = ICOS_beam(@ICOS_Model6, P);
    IB.Sample('beam_samples', Opt.Nsamples, ...
      'ICOS_passes', Opt.ICOS_passes, 'opt_n', opt_n, ...
      'n_optics', n_optics, 'Track_Power', Track_Power, ...
      'Herriott_passes', 1, ...
      'mnc', IBmnc, IBopt{:});
    if opt_n == n_optics
      ff2 = IB.Integrate;
      if ~isempty(ff)
        delete(ff);
        drawnow;
      end
      ff = ff2;
    end
    IB.savefile;
    % save(ofile, 'IB');
    % fprintf(1, 'ICOS_beam %d Saved result to %s\n', i, ofile);
    Pwr = IB.PowerSummary;
    IS.res2(i).NH = Pwr.NH;
    if isfield(Pwr,'max_pwr')
      IS.res2(i).max_pwr = Pwr.max_pwr;
    elseif ~isfield(IS.res2,'max_pwr')
      IS.res2(i).max_pwr = [];
    end
    IS.savefile;
  else
    fprintf(1, 'Skipping result %d\n', i);
  end
