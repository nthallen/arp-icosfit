%%
% This script needs to run after scratch_151104.m
% It selects promising solutions and sends them through ICOS_search
% to get complete configurations. From there, we could analyze.
RLmin = [Summary.RLmin];
dL = [Summary.dL];
ddL = diff(dL);
%%
figure;
plot(RLmin,dL,'.');
ylim([-1.5 1.5]);
%%
plot(RLmin,ddL,'.');
ylim([0 2]);
xlim(min(RLmin)+[-0.2 5]);
%%
ri = find(RLmin < 20 & ddL < 1.1);
Sums = Summary(ri);
%%
Ifocuses = zeros(1,length(Sums));
for i = 1:length(Sums)
  %%
  SP.R1 = Sums(i).R1;
  SP.R2 = Sums(i).R2;
  SP.L = L;
  SP.Rw1 = 1.1*B/2;
  phi = pi/Sums(i).m;
  Phi = Sums(i).k*phi;
  r_max = sqrt(rdtanth*sqrt(((SP.R2-L).*SP.R1.^2*L)./((SP.R1-L).*(SP.R1+SP.R2-L))));
  if isnan(r_max)
    error('r_max %d is NaN', i);
  end
  r_max = min(r_limit,r_max);
  r2 = abs(r_max.*cos(Phi).*SP.R2./(SP.R2-L));
  w2 = abs(r_max*sin(Phi)*cos(Phi).*SP.R2./(SP.R2-L));
  s1 = w2/L;
  SP.RL = SP.Rw1/s1;
  Res = exparam(SP);
  check_params(i, Res);
  mnc = sprintf('sropt_g.L%d.%d',L,i);
  IS = ICOS_search('mnc', mnc,'R1',SP.R1,'R2',SP.R2,'L',L,'RR1',Res.RR1,'Rw1',SP.Rw1);
  %%
  IS.search_ICOS_RIM;
  %%
  IS.search_focus2('max_lenses',2);
  Ifocuses(i) = length(IS.res2);
end
%%
for i=1:length(Sums)
  fname = sprintf('IS_sropt_g.L%d.%d.mat',L,i);
  load(fname);
  IS.search_focus2('max_lenses',2);
end
%%
% This displays all the solutions
for i=1:length(Sums)
  fname = sprintf('IS_sropt_g.L%d.%d.mat',L,i);
  load(fname);
  for j=1:length(IS.res2)
    PM = ICOS_Model6(render_model(IS.res2(j),'ICOS_passes_per_injection', 1,'view',[0 90]));
    title(sprintf('%s solution %d', strrep(fname,'_','\_'), j));
    pause;
  end
end
%%
% Collect the results
res = collect_results('files','IS_sropt_g.*.mat','exclude','NH', ...
           'exclude','max_pwr','Ltot',[]);
[~,index] = sortrows([res.Ltot].'); res = res(index); clear index
%%
% This approach to parallel processing is flawed. Can't have more than
% one thread working on one IS object at a time. However, it can be
% cleaned up after the fact.
ncores = 4;
core = 1;
for i=core:ncores:length(res)
  fname = sprintf('IS_%s.mat', res(i).mnc);
  IS = load(fname);
  IS = IS.IS;
  IS.analyze('select', res(i).index);
end
%%
% Oops, my parallel processing strategy was flawed. Need to fix up
% files = dir('IS_sropt_g.L50.*.mat');
% files = { files.name };
% %%
% isIB = regexp(files,'(\dx\d)|(_save)');
% cellisempty = @(x) isempty(x{1});
% isIB2 = arrayfun(cellisempty, isIB);
% files = files(isIB2)';
% %%
% for i=1:length(files)
%   %%
%   load(files{i}); % defines IS
%   %%
%   for j=1:length(IS.res2)
%     if isempty(IS.res2(j).NH)
%       IBfnamepat = sprintf('IS_%s.%d_*.mat', IS.ISopt.mnc, IS.res2(j).Nres2);
%       IBfile = dir(IBfnamepat);
%       if length(IBfile) ~= 1
%         fprintf(1,'Pattern "%s" yielded %d results\n', IBfnamepat, ...
%           length(IBfile));
%       else
%         load(IBfile.name);
%         Pwr = IB.PowerSummary;
%         IS.res2(j).NH = Pwr.NH;
%         if isfield(Pwr,'max_pwr')
%           IS.res2(j).max_pwr = Pwr.max_pwr;
%         elseif ~isfield(IS.res2,'max_pwr')
%           IS.res2(j).max_pwr = [];
%         end
%       end
%     end
%   end
%   IS.savefile;
% end
%%
% Collect post-analysis results
res = collect_results('files','IS_sropt_g.*.mat','Rr1',[],'R1',[],'D1',[],'Ltot',[]);
%[~,index] = sortrows([res.Ltot].'); res = res(index); clear index
%%
scratch_plot(res, 'H_loss','max_pwr');
%%
scratch_plot(res, 'I_loss','max_pwr');
%%
scratch_plot(res, 'D_loss','max_pwr');
%%
% Sanity check on power losses
for i=1:length(res)
  res(i).cum_trans = (1-res(i).H_loss/100) * ...
    (1-res(i).I_loss/100) * ...
    prod(1 - res(i).F_loss/100) * ...
    (1-res(i).D_loss/100);
  res(i).proj_pwr = res(i).NH * res(i).cum_trans;
end
%%
scratch_plot(res,'cum_trans','max_pwr');
%%
scratch_plot(res,'proj_pwr','max_pwr');
%%
V = polyfit([res.proj_pwr],[res.max_pwr],1);
%%
% Sanity check that NH ~= pi/asin(Rw1/Rr1);
Rw1 = [res.Rw1];
Rr1 = [res.Rr1];
NHphi = floor(pi./asin(Rw1./Rr1));
NH = [res.NH];
figure; plot(NHphi,NH-NHphi,'.');
xlabel('NH \phi');
ylabel('NH');
%%
% Check if H_loss is related to radial margin
H_r = [res.RD1]*2.54/2;
Hmargin = H_r - [res.Rr1];
figure; plot(Hmargin,[res.max_pwr],'.');
xlabel('Margin, cm');
ylabel('max\_pwr');
%%
% Check if H_loss is related to radial margin on ICOS mirror:
I_r = [res.D1]*2.54/2;
HImargin = I_r - [res.r1];
figure; plot(HImargin,[res.max_pwr],'.');
xlabel('HIMargin, cm');
ylabel('max\_pwr');
