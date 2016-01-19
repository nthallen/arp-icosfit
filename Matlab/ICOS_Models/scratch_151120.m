%% scratch_151120
% Developing ICOS_sr_search
SR = ICOS_sr_search('mnc','sr.L50t','B',0.4,'Rw1',0.2,'mk',[31,2],'SolScale', 0.5);
SR.enumerate;
%%
SR.design_tolerance;
SR.build_tolerance;
%%
SR.explore_build_tolerance;
%%
SR.explore_build_tolerance(2);
%%
SR.explore_build_tolerance(3);
%%
SR.explore_build_tolerance(4);
%%
SR.savefile;
%%
SR.focus('focus_visible',0);
%% Try the same parameters, but push Rw1 wider to see if we can cut down
% on H loss
SR = ICOS_sr_search('mnc','sr.L50.Rw1.5','B',0.4,'Rw1',0.6);
SR.enumerate;
SR.design_tolerance;
SR.build_tolerance;
%%
SR.explore_build_tolerance;
%%
SR.explore_build_tolerance(2);
%%
SR.savefile;
SR.focus('focus_visible', 0);
%%
files = dir('IS_sr.L50a.*.mat');
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
RLmin = [SR.Summary.RLmin];
dL = [SR.Summary.dL];
ddL = diff(dL);
dBL = [SR.Summary.dBL];
ddBL = diff(dBL);
%%
figure;
plot(ddL,ddBL,'.');
xlabel('Design \Delta L cm');
ylabel('Build \Delta L cm');
%%
figure;
scatter(ddL,ddBL,[],RLmin);
xlabel('Design \Delta L cm');
ylabel('Build \Delta L cm');
colorbar
%%
figure;
v = ddL < 2000;
scatter(RLmin(v),ddBL(v),[],ddL(v));
xlabel('RLmin cm');
ylabel('Build \Delta L cm');
colorbar
%%
figure;
v = ddL < 2000;
m = [SR.Summary.m];
scatter(RLmin(v),ddBL(v),[],m(v));
xlabel('RLmin cm');
ylabel('Build \Delta L cm');
cb = colorbar;
ylabel(cb, 'm');
%%
v = ddL < 2000;
k = [SR.Summary.k];
vi = find(v);
%%
figure;
plot(RLmin(v),ddBL(v),'.');
%set(gca,'UserData',xtra_text);
xlabel('RLmin cm');
ylabel('Build \Delta L cm');
hdt = datacursormode;
set(hdt,'UpdateFcn',{@data_cursor_text_func,RLmin,ddBL,m,k})
%%
figure;
x = rand(30,1);
y = rand(30,1);
v = x > .5 & y > .5;
h = plot(x,y,'.',x(v),y(v),'o');
arg3 = 1;
set(h, 'buttondownfcn', @(src,ev) test_button_down(src,ev,3));
arg3 = 2;
tbd2 = @(src,ev) test_button_down(src,ev,arg3);
set(gca, 'buttondownfcn', tbd2);
%% ICOS_sr_search exploring alternatives on 31,2
for i=1:10
  mnc = sprintf('L50b_31.2.%d', i);
  SR = ICOS_sr_search('mnc',mnc,'B',0.4,'Rw1',0.2,'mk',[31,2],'SolScale', i/10);
  SR.enumerate;
  %%
  SR.design_tolerance;
  SR.build_tolerance;
  SR.savefile;
  SR.focus('focus_visible',0);
end
%%
SR.explore_build_tolerance;
%%
SR.explore_build_tolerance(2);
%%
SR.explore_build_tolerance(3);
%%
SR.explore_build_tolerance(4);
%%
SR.savefile;
%%
SR.focus('focus_visible',0);
%%
files = dir('IS_L50b_31.2.*.mat');
for i=1:10
  load(files(i).name);
  fprintf(1,'Analyzing %s\n', IS.ISopt.mnc);
  IS.analyze('Nsamples',400);
  close all
end
fprintf(1,'Done\n');
%%
res = collect_results('files','IS_L50b*.mat');
%%
for i = 1:10
  mnc = sprintf('L50b_31.2.%d.1_31.2',i);
  ni = reshape(find(strcmp({res.mnc},mnc)),1,[]);
  for nii=ni
    res(nii).SolScale = i/10;
  end
end
%%
plot([res.SolScale],[res.max_pwr],'.');
xlabel('SolScale');
ylabel('max\_pwr');
%%
plot([res.SolScale],[res.RL],'.');
xlabel('SolScale');
ylabel('RL');
%%
SolScale = [res.SolScale];
max_pwr = [res.max_pwr];
SSbins = unique(SolScale);
max_pwr_bins = zeros(size(SSbins));
for i=1:length(SSbins)
  v = SolScale == SSbins(i);
  max_pwr_mean(i) = mean(max_pwr(v));
  max_pwr_std(i) = std(max_pwr(v));
end
plot(SolScale, max_pwr, '.', SSbins, max_pwr_mean, '*-');
xlabel('SolScale parameter');
ylabel('Max Power');
xlim([0.05 1.05]);
%%
load('IS_L50b_31.2.4.1_31.2.mat');
v = find([IS.res2.sel]);
LS = zeros(length(v),4);
for i=1:length(v)
  LS(i,:) = [IS.res2(v(i)).Lens_Space, IS.res2(v(i)).detector_spacing, ...
    IS.res2(v(i)).max_pwr];
end

