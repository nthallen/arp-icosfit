function print_figure_v1(ax, fname, target_dim, t)
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
  
  fig = get(ax,'parent');
  set(fig,'InvertHardcopy','on');
  set(fig,'PaperUnits', 'inches');
  papersize = get(fig, 'PaperSize');
  % Center the image on the page:
  left = (papersize(1)- width)/2;
  bottom = (papersize(2)- height)/2;
  myfiguresize = [left, bottom, width, height];
  set(fig,'PaperPosition', myfiguresize);
  
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
  drawnow;
  
  set(ax,'xlimmode','manual','ylimmode','manual','zlimmode','manual');
  
  print(fig, '-dpng','-r300',fname);
  
  for i=1:length(t)
    set(t(i),'fontsize', tpixels(i));
  end
  % set(h12,'linewidth',h12w);
  % set(h23,'linewidth',h23w);
  drawnow;
  
