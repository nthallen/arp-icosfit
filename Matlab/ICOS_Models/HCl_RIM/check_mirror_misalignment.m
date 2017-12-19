%%
% Check mirror misalignment
cd C:\Users\nort.ARP\Documents\SW\arp-icosfit\Matlab\ICOS_Models\HCl_RIM
load('IB_HCl_RIM100_B2_Rw3_rd08_th14.0.14_200_1000_100_40.7.1_50x100.mat');
%%
P = IB.P;
P.y0 = 0; % On-axis alignment
P.dy = 0; 
P.dz = 0;
P.HR = 0; % No Herriott reflection
P.M2dy = 0;
%%
% Pick a mirror deflection (M2dy)
% See how it affects power to the detector, then try to compensate
% by adjusting y0 and/or dy
P.visible = 1;
P.evaluate_endpoints = 6;
PM = ICOS_Model6(P);
%%
P.D_l = 0.1;
max_non_parallel = .01*2.54;
max_M2dy = max_non_parallel/(2*P.r2);
P.visible = 0;
P.plot_endpoints = 0;
N = 100;
PM = ICOS_Model6(P,'M2dy',(0:N)*max_M2dy/N);
PM.plot_results('total_power');
%%
P.M2dy = max_M2dy;
P.plot_endpoints = 0;
P.visible = 1;
PM = ICOS_Model6(P);
%%
P.M2dy = max_M2dy;
P.plot_endpoints = 0;
P.visible = 0;
dy = 0.8;
ddy = 0.02;
PM = ICOS_Model6(P, 'y0', [-5:5]*dy/5, 'dy', [-5:5]*ddy/5);
PM.plot_results('total_power');
shg;
%%
P.y0 = 0.7;
P.dy = .007;
P.visible = 1;
PM = ICOS_Model6(P);
P.y0 = 0; 
P.dy = 0;
