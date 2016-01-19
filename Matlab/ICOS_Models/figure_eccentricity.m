function figure_eccentricity
%% This code is extracted from ICOS_search.search_ICOS_RIM
load('IS_L50c.4.mat');
P = render_model(IS.res2(12));
P.stop_ICOS = 0;
P.visible = 0;
P.visibility = 0;
P.focus = 0;
P.HR = 0;
P.ICOS_passes_per_injection = 100;
P.max_rays = 3000;
P.plot_endpoints = 0;
P.evaluate_endpoints = 3;
P.skip.overlap = 1;
P.skip.total_power = 1;
P.skip.mean_angle = 1;
P.skip.RIM_passes = 1;
delta = 2.5e-3;
best_eccentricity = 2;
eccentricity = 1;
iteration = 1;
while delta > 1e-4 && eccentricity < best_eccentricity && iteration < 10
  best_eccentricity = eccentricity;
  PM = ICOS_Model6(P,'dy',P.dy+linspace(-delta,delta,21),'dz',P.dz+linspace(-delta,delta,21));
  PM = PM.clean_results;
  [new_dy,new_dz,eccentricity,ri,ci] = PM.identify_minimum('eccentricity');
  if isempty(new_dy) || isempty(new_dz)
    delta = delta*10;
  elseif eccentricity < best_eccentricity
    P.dy = new_dy;
    P.dz = new_dz;
    if ri ~= 1 && ri ~= size(PM.Results.eccentricity,1) && ...
        ci ~= 1 && ci ~= size(PM.Results.eccentricity,2)
      delta = delta/2;
    end
  end
  iteration = iteration + 1;
end
%%
figure;
PM.plot_results('eccentricity');
view([-77.5, 18.8]);
set(gca,'fontsize',12,'fontweight','bold');
%%
 
print_figure(gca, 'figure_eccentricity.png', [3 3],gca);

function print_figure(ax, fname, target_dim, t)
  target_width=target_dim(1);
  target_height=target_dim(2);
  
  target_aspect_ratio = target_height/target_width;
  set(ax,'units','inches');
  P = get(ax,'position');
  source_aspect_ratio = P(4)/P(3);
  if source_aspect_ratio >= target_aspect_ratio
    % fit to target height
    height = target_height;
    width = height/source_aspect_ratio;
    print_scale = height/P(4);
  else
    width = target_width;
    height = width * source_aspect_ratio;
    print_scale = width/P(3);
  end
  
  set(gcf,'InvertHardcopy','on');
  set(gcf,'PaperUnits', 'inches');
  papersize = get(gcf, 'PaperSize');
  % Center the image on the page:
  left = (papersize(1)- width)/2;
  bottom = (papersize(2)- height)/2;
  myfiguresize = [left, bottom, width, height];
  set(gcf,'PaperPosition', myfiguresize);
  
  % Now change the font size to map the size we will be scaled to
  tpixels = zeros(size(t));
  for i=1:length(t)
    set(t(i),'fontunits','pixels');
    tpixels(i) = get(t(i),'fontsize');
    set(t(i),'fontsize', tpixels(i)*print_scale);
  end
  
%   h12w = get(h12,'linewidth');
%   set(h12,'linewidth',h12w*print_scale);
%   
%   h23w = get(h23,'linewidth');
%   set(h23,'linewidth',h23w*print_scale);
  
  set(ax,'xlimmode','manual','ylimmode','manual','zlimmode','manual');
  
  print(fname,'-dpng','-r300');
  
  for i=1:length(t)
    set(t(i),'fontsize', tpixels(i));
  end
  % set(h12,'linewidth',h12w);
  % set(h23,'linewidth',h23w);
