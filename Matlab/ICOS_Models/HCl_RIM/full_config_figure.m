%%
load('IB_HCl_RIM100_B2_Rw3_rd08_th14.0.14_200_1000_100_40.7.2_50x100.mat')
IB.P.r1 = 2.54;
IB.P.r2 = 2.54;
IB.P.Hr = 1.5*2.54;
IB.P.ICOS_passes_per_injection = 1;
IB.draw;
%%
padding = [-2 2];
xy = 1.6*2.54*[-1 1];
set(gca,'zlim',xy,'ylim',xy+padding,'xlim',[-32 56]+padding);
set(gca,'position',[0 0 1 1]);
view(0,90)
%%
set(gca,'visible','off');
set(gcf,'color','white');
%%
print(gcf, 'full_config_figure.png', '-dpng');
%%
IB.P.visibility = [0 0 0];
IB.draw;
view(0,90);
