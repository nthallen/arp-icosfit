function run_151212(core, cores)
files = dir('IS_L50b_31.2.*.mat');
for i=core:cores:length(files)
  load(files(i).name);
  fprintf(1,'Analyzing %s\n', IS.ISopt.mnc);
  IS.analyze('Nsamples',400);
  close all
end
fprintf(1,'Done\n');
