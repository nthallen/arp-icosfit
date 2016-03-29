function rt_quick_plot(res, fld1, fld2)
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
  end
