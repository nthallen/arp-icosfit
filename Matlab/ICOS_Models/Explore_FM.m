%% Explore Focus_Model
cd C:\Data\ES96\ICOS\3D
P = Focus_Model.props;

PM = Focus_Model(P);
%%
P.visible = true;
P.pause = false;
P.plot_endpoints = 0;
PM = Focus_Model(P,'D_X',P.L1_EFL * linspace(1,1.2,10));
%%
PM.plot_results('max_radius');
%%
[X, ~, max_radius] = PM.identify_minimum('max_radius');
P.plot_endpoints = 0;
P.visible = false;
PM = Focus_Model(P,'D_X',X + 0.2 * P.L1_EFL * linspace(-1,1,20));
PM.plot_results('max_radius');
[X, ~, max_radius] = PM.identify_minimum('max_radius');

%%
P.D_X = X;
P.dz = 0;
P.visible = true;
PM = Focus_Model(P);
%%
S.Parallel_FL = X;
fprintf(1,'Parallel focal point is %.2f or %.2f past EFL\n', X, X-P.L1_EFL);
%%
% Now let's look at what divergence does to our focus
P.dz = 0.036; % This number was taken from Explore_IM4
PM = Focus_Model(P);

%%
P.visible = false;
P.pause = true;
P.plot_endpoints = 2;
PM = Focus_Model(P,'D_X',P.D_X * linspace(1,1.2,10));
%%
PM.plot_results('max_radius');
%%
[X, ~, max_radius] = PM.identify_minimum('max_radius');
P.plot_endpoints = 2;
PM = Focus_Model(P,'D_X',X + 0.2 * P.L1_EFL * linspace(-1,1,20));
PM.plot_results('max_radius');
[X, ~, max_radius] = PM.identify_minimum('max_radius');
%%
P.visible = true;
P.pause = false;
P.D_X = X;
PM = Focus_Model(P);
%%
S.Divergent_FL = X;
fprintf(1,'Divergent focal point is %.2f or %.2f past EFL\n', X, X-P.L1_EFL);
%%
PM.M.plot_endpoints(2);

%% Now let's look at skew

