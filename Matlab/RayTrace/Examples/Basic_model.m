%% Basic example of opt_model use:
n_ZnSe = 2.4361;
n_air = 1;
P.max_rays = 100;
P.r = 3*2.54/2;
P.RC = 150;
P.CT = 0.2; % mirror center thickness
P.d = 50; % mirror spacing
P.T = 250e-6; % transmittance
P.y0 = 2.9;
P.dy = .04;
P.dz = .03;
P.vis = true;

M = opt_model(2, P.max_rays);
M.visible = P.vis;
M.Optic{1} = HRmirror('M1', P.r, P.RC, P.CT, P.T, 1-P.T, [0 0 0], [1 0 0], n_ZnSe, n_air, P.vis);
M.Optic{2} = HRmirror('M2', P.r, P.RC, P.CT, P.T, 1-P.T, [P.d 0 0], [-1 0 0], n_ZnSe, n_air, P.vis);
M.push_ray(opt_ray([-1, P.y0, 0], [1, -P.dy, P.dz]), 0, 0, 1);
M.propagate;
