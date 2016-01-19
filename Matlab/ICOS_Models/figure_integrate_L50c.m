%% figure_integrate_L50c.m
load('IB_L50c.4.14_50x100.mat');
%%
ff = IB.Integrate;
ax = ones(size(ff));
for i=1:length(ff)
  ax(i) = findobj(ff(i),'type','axes');
  set(ax(i),'fontsize',12,'fontweight','bold');
end
title(ax(1),'High Resolution');
title(ax(2),'Detector Resolution');
%%
for i=1:length(ff)
  print_figure_v1(ax(i),sprintf('figure_integrate_L50c_%d.png', i),[3 3], ax(i));
end
