%%
P = ICOS_Model6.props;
LT = P.LensTypes;
Lnames = fieldnames(LT);
%%
nL = length(Lnames);
EFL = zeros(nL,1);
r = zeros(nL,1);
R1 = zeros(nL,1);
for i=1:nL
  Lens = LT.(Lnames{i});
  fl = Lens.EFL;
  if strcmp(Lens.type,'positive_meniscus')
    if fl < 0
      fprintf(1, 'Lens %s is PM with negative EFL\n', Lnames{i}, fl);
      fl = -fl;
    end
  elseif strcmp(Lens.type,'negative_meniscus')
    if fl > 0
      fprintf(1, 'Lens %s is NM with posative EFL\n', Lnames{i}, fl);
      fl = -fl;
    end
  else
    fprintf(1,'Lens %s has unrecognized type %s\n', Lnames{i}, Lens.type);
  end
  EFL(i) = fl;
  r(i) = Lens.r;
  R1(i) = Lens.R1;
end
%%
% Create a new postive meniscus lens with radius 1.25", say, with
% EFL of 2:5 inches. Use the thick lens formulat from the ISP Optics
% catalog. Start with a range of CT values, then interpolate to match
% the target ET value.
r = 1.25 * 2.54;
ET = 0.4;
M = 0.9;
n = 2.4361;
for EFLin = [2,3,4,5]
  EFL = EFLin*2.54;
  name = sprintf('P.LensTypes.Cust_ZC_PM_%d_%d', floor(r*20), floor(EFL*10));
  CTv = ET+[0:.01:2];
  R1 = M*EFL;
  R2v = (n-1)*EFL*((1-((CTv*(n-1))/(M*EFL*n)))/(((n-1)/M)-1));
  sag1 = R1 - sqrt(R1.^2-r^2);
  sag2 = R2v - sqrt(R2v.^2-r^2);
  ETv = CTv-sag1+sag2;
  CT = interp1(ETv,CTv,ET);
  R2 = interp1(ETv,R2v,ET);
  fprintf(1,'%s.type = ''positive_meniscus'';\n', name);
  fprintf(1,'%s.r = %.4f;\n', name, r);
  fprintf(1,'%s.R1 = %.4f;\n', name, R1);
  fprintf(1,'%s.R2 = %.4f;\n', name, R2);
  fprintf(1,'%s.CT = %.4f;\n', name, CT);
  fprintf(1,'%s.EFL = %.4f;\n', name, EFL);
end  
