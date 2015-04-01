%%
P = ICOS_Model3.props;
P.visible = true;
P.plot_endpoints = 0;
P.pause = 1;
P.title = 'Scan';
P.fmt_L2_X = 'L2\\_X = %.1f cm';
P.fmt_L1_R2 = 'L1\\_R2 = %.2f cm';
% P.avifile = 'L2_sweep.avi';
P.dz = 0.01;
P.max_rays = 200;
PM = ICOS_Model3(P,'L2_X',linspace(57,62,20));
%PM = ICOS_Model3(P,'L1_R2',linspace(27.9,29.9,10));
%PM = ICOS_Model3(P);
%%
PM.M.plot_endpoints(5);
%%
M = P.L1_R1/P.L1_EFL;
n_ZnSe = 2.4361;
n_air = 1;
n_int = n_ZnSe;
EFL = P.L1_EFL;
R1 = P.L1_R1;
R2 = (n_int-1)*EFL* ...
        (1-(P.L1_CT*(n_int-1)/(M*EFL*n_int)))/...
        ((n_int-1)/M - 1);
%%
R2a = (EFL * (n_int-1))/((n_int-1)/M-1);