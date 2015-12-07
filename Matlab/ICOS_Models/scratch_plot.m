function scratch_plot(res, fld1, fld2)
  figure;
  plot([res.(fld1)], [res.(fld2)], '.');
  xlabel(strrep(fld1,'_','\_'));
  ylabel(strrep(fld2,'_','\_'));
  
