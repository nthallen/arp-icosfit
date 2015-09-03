%%
% Analyze relationship between number of samples and standard deviation
%
% Tentative conclusion: Taking many samples and analyzing does not appear
% to estimate properly. Combining results from smaller sample sizes
% seems to behave normally. There seems to be a problem with how the larger
% samples are analyzed that somehow suppresses variation. My guess is
% roundoff error when summing the power on the detector.
mnc = 'R28_w25_r24d';
filename = sprintf('IS_%s.mat', mnc);
load(filename);
Nres2 = [IS.res2.Nres2];
filepat = sprintf('IS_%s.*_*x*.mat', mnc);
files = dir(filepat);

Nsamples = zeros(length(files),1);
Nres2 = zeros(length(files),1);
Npasses = zeros(length(files),1);

for i=1:length(files)
  pat = sprintf('IS_%s.%%d_%%dx%%d.mat', mnc);
  res = sscanf(files(i).name, pat);
  Nres2(i) = res(1);
  Npasses(i) = res(2);
  Nsamples(i) = res(3);
end
%%
clear Power;
for i=length(files):-1:1
  load(files(i).name);
  Pwr = IB.PowerSummary;
  Pwr.Nres2 = Nres2(i);
  Pwr.Nsamples = Nsamples(i);
  Pwr.Npasses = Npasses(i);
  Power(i) = Pwr;
end
clear IB
%% Look at the mean and std of power
Nsamples = [Power.Nsamples];
NS = unique(Nsamples);
max_pwr = [Power.max_pwr];
for ns = NS
  v = Nsamples == ns;
  Pmean = mean(max_pwr(v));
  Pstd = std(max_pwr(v));
  fprintf(1,'NS=%d mean Power: %f std Power: %f  std mean: %f\n', ...
    ns, Pmean, Pstd, Pstd/sqrt(sum(v)));
end
%% Now let's combine the NS(1) runs in pairs of two, recalculate max_pwr
%  and see how that affects mean and std
ns = NS(1);
vi = find(Nsamples == ns);
v2 = find(Nsamples == NS(2));
std1 = std(max_pwr(vi));
std2 = std(max_pwr(v2));
fprintf(1, 'std(NS=%d) = %f\n', NS(1), std1);
fprintf(1, 'std(NS=%d) = %f\n', NS(2), std2);
%%
N_permutations = 100;
std1c = zeros(N_permutations,1);
for p_i = 1:N_permutations
  %%
  vi = vi(randperm(length(vi)));
  npairs = floor(length(vi)/2);
  mean_pwr = zeros(npairs,1);
  combined_pwr = zeros(npairs,1);
  for i=1:npairs
    %%
    load(files(vi(2*i-1)).name);
    IB1 = IB;
    load(files(vi(2*i)).name);
    Int = [IB1.Res.Int(2); IB.Res.Int(2)]; % Implicitly taking Dyz == 0.2
    %%
    Yr = [min(Int(1).Yr(1),Int(2).Yr(1)) max(Int(1).Yr(2),Int(2).Yr(2))];
    Zr = [min(Int(1).Zr(1),Int(2).Zr(1)) max(Int(1).Zr(2),Int(2).Zr(2))];
    %%
    acc = zeros(Yr(2)-Yr(1)+1,Zr(2)-Zr(1)+1);
    % image(Y,Z,Int.img') puts Y on the X-axis and Z on Y-axis
    % Int.img(yi,zi) is how it is addressed.
    %%
    for j=1:2
      %%
      Yi = (Int(j).Yr(1):Int(j).Yr(2)) - Yr(1)+1;
      Zi = (Int(j).Zr(1):Int(j).Zr(2)) - Zr(1)+1;
      acc(Yi,Zi) = acc(Yi,Zi) + Int(j).img;
    end
    acc = acc/2;
    mean_pwr(i) = (max_pwr(vi(2*i-1)) + max_pwr(vi(2*i)))/2;
    combined_pwr(i) = max(max(acc));
  end
  %
  % figure;
  % plot(1:npairs,mean_pwr-combined_pwr,'+');
  % legend('mean pwr');
  %
  std1c(p_i) = std(combined_pwr);
  fprintf(1, 'std(N=%d,2): %f\n', NS(1), std1c(p_i));
end
%%
% Look at the 200-sample runs
% Divide each into two 100-sample runs and re-Integrate
v200 = find(Nsamples == 200);
Nv200 = length(v200);
Psum(Nv200) = struct('Pmax200',[],'Pmax100',[]);
for v200i=1:Nv200;
  load(files(v200(v200i)).name);
  Pwr = IB.PowerSummary;
  Psum(v200i).Pmax200 = Pwr.max_pwr;
  Pwrs = [0 0];
  for j=[0 1]
    load(files(v200(v200i)).name);
    resi = (1:length(IB.Res.Sample))' <= IB.Res.N;
    if j
      vs = resi & (IB.Res.Sample <= 100);
    else
      vs = resi & (IB.Res.Sample > 100);
    end
    IB.IBP.beam_samples = 100;
    IB.Res.D = IB.Res.D(vs,:);
    IB.Res.E = IB.Res.E(vs,:);
    IB.Res.P = IB.Res.P(vs,:);
    IB.Res.NPass = IB.Res.NPass(vs,:);
    IB.Res.Sample = IB.Res.Sample(vs)-j*100;
    IB.Res.NPasses = IB.Res.NPasses((1:100)+j*100);
    IB.Res.N = sum(vs);
    IB.Res.Int(end).img = []; % Just integrate for the big detector
    ff = IB.Integrate;
    title(sprintf('%d/%d split %d', v200i, length(v200), j));
    delete(ff(1));
    Pwr = IB.PowerSummary;
    Pwrs(j+1) = Pwr.max_pwr;
  end
  Psum(v200i).Pmax100 = Pwrs;
  fprintf(1, '%d: %f vs %f %f\n', v200i, Psum(v200i).Pmax200, ...
    Psum(v200i).Pmax100(1), Psum(v200i).Pmax100(2));
end
