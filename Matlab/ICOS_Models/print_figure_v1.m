function print_figure_v1(ax, fname, target_dim, t)
  % print_figure_v1(ax, fname, target_dim, t);
  % ax: axes or figure to print
  % fname: name of output png file
  % dim: 1x2, [w,h] in inches
  % t: vector of text objects to be scaled for printing
  target_width=target_dim(1);
  target_height=target_dim(2);
  
  target_aspect_ratio = target_height/target_width;
  if strcmp(get(ax,'type'),'figure')
    fig = ax;
  else
    fig = get(ax,'parent');
  end
  set(fig,'units','inches');
  P = get(fig,'position');
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
    if ~strcmp(get(t(i),'type'),'colorbar')
      set(t(i),'fontunits','pixels');
    end
    tpixels(i) = get(t(i),'fontsize');
    set(t(i),'fontsize', tpixels(i)*print_scale);
  end
  
  %   h12w = get(h12,'linewidth');
  %   set(h12,'linewidth',h12w*print_scale);
  %
  %   h23w = get(h23,'linewidth');
  %   set(h23,'linewidth',h23w*print_scale);
  drawnow;
  
  axs = findobj(fig,'type','axes');
  set(axs,'xlimmode','manual','ylimmode','manual','zlimmode','manual');
  
  print(fig, '-dpng','-r300',fname);
  
  for i=1:length(t)
    set(t(i),'fontsize', tpixels(i));
  end
  % set(h12,'linewidth',h12w);
  % set(h23,'linewidth',h23w);
  drawnow;
  
