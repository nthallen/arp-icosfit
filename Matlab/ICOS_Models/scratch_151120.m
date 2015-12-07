%% scratch_151120
% Developing ICOS_sr_search
SR = ICOS_sr_search('mnc','sr.L50');
SR.enumerate;
SR.design_tolerance;
SR.build_tolerance;
%%
SR.explore_build_tolerance;
%%
SR.explore_build_tolerance(2);
%%
SR.savefile;
SR.focus;
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
