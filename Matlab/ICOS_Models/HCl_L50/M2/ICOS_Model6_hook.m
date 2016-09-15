function M = ICOS_Model6_hook( M, ~ )
  % M = ICOS_Model6_hook(M,P);
  %   Used here to tweak the model M before running analyses
  M.Optic{2}.T = 1;
  M.Optic{2}.R = 0;
  M.Optic{2}.n_int = 1;
  M.Optic{2}.Surface{1}.T = 1;
  M.Optic{2}.Surface{1}.R = 0;
  M.Optic{2}.Surface{1}.n_int = 1;
  M.Optic{2}.Surface{2}.T = 1;
  M.Optic{2}.Surface{2}.R = 0;
  M.Optic{2}.Surface{2}.n_int = 1;
end

