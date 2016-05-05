function rt_quick_plot(res, fld1, fld2, fld3)
  % rt_quick_plot(res, fld1, fld2);
  % rt_quick_plot(res, fld1, fld2, fld3);
  % Examples:
  %   rt_quick_plot(res, 'R1', 'L');
  %   rt_quick_plot(res, 'R1', 'R2', 'L');
  figure;
  Xdata = [res.(fld1)];
  if nargin == 2
    plot(Xdata,'.');
    ylabel(strrep(fld1,'_','\_'));
  elseif nargin == 3
    Ydata = [res.(fld2)];
    plot(Xdata,Ydata,'.');
    xlabel(strrep(fld1,'_','\_'));
    ylabel(strrep(fld2,'_','\_'));
  elseif nargin == 4
    Ydata = [res.(fld2)];
    Cdata = [res.(fld3)];
    scatter(Xdata,Ydata,[],Cdata);
    xlabel(strrep(fld1,'_','\_'));
    ylabel(strrep(fld2,'_','\_'));
    cb = colorbar;
    cb.Label.String = fld3;
  end
