function feedback_analysis(varargin)
% feedback_analysis([scans,] [refbase ,] base);
% 
% feedback_analysis takes one to three arguments:
%   An optional list of scan numbers to use in the analysis
%   An optional reference base name of an icosfit for comparison
%   The base name of an icosfit run to fully analyze.
base = '';
refbase = '';
scans = [];
switch length(varargin)
  case 1,
    base = varargin{1};
  case 2,
    if isnumeric(varargin{1})
      scans = varargin{1};
      base = varargin{2};
    else
      refbase = varargin{1};
      base = varargin{2};
    end
  case 3,
    scans = varargin{1};
    refbase = varargin{2};
    base = varargin{3};
  otherwise
    help feedback_analysis
    return;
end
[line_pos,S] = get_line_position(base);
if ~isempty(refbase)
  [line_pos0,S0] = get_line_position(refbase);
  osc0 = fit_comp([], S0, line_pos0, scans);
  osc = fit_comp(S0, S, line_pos, scans);
  figure;
  lnums = 1:length(osc);
  plot(lnums,osc0,'rs-',lnums,osc,'*-b');
  set(gca,'YGrid','On','XDir','Reverse','xlim',[0.75 lnums(end)+0.25]);
  title(sprintf('%s/%s: Residual Oscillation by Line',getrun,base));
  legend('Standard','Feedback');
else
  osc = fit_comp([], S, line_pos, scans);
  figure;
  lnums = 1:length(osc);
  plot(lnums,osc,'rs-');
  set(gca,'YGrid','On','XDir','Reverse','xlim',[0.75 lnums(end)+0.25]);
  title(sprintf('%s/%s: Residual Oscillation by Line',getrun,refbase));
end


function osc_out = fit_comp(S0, S, line_pos, scans)
if isempty(scans)
  scans = S.scannum;
end
[scans,vs] = scan_intersect(scans, S.scannum);
PTE = load(S.PTEfile);
[Pscans,vp] = scan_intersect(scans, PTE(:,1));
if length(Pscans) ~= length(scans)
  error('Scans present in S are not present in PTE file');
end
PTE = PTE(vp,:);
scan10 = fastavg(scans,10);
if ~isempty(S0)
  [scans0,vs0] = scan_intersect(scans, S0.scannum);
  scan010 = fastavg(scans0,10);
end

x = (PTE(1,4):max(S.SignalRegion(:,2)))';
X = (x' - PTE(1,4))/1000;
row = ones(1,length(x));
col = ones(length(scans),1);
a = PTE(:,7);
b = PTE(:,6);
c = PTE(:,5);
d = PTE(:,8);
tau = PTE(:,9);
d2 = PTE(:,10);
tau2 = PTE(:,11);
fn = a*(X.^2) + b*X + c*row + (d*row).*exp(-(1./tau)*X) + ...
    (d2*row).*exp(-(1./tau2)*X);
[XX,YY] = meshgrid(x,scans);

osc = zeros(1,size(S.Chi,2));
for loi=1:size(S.Chi,2)
  f = figure;
  p = get(f,'Position');
  new_h = 500;
  new_w = 1500;
  p = [ p(1)-(new_w-p(3))/2 p(2)-(new_h-p(4)) new_w new_h];
  set(f,'Position',p)
  xsp = .08;
  if isempty(S0)
    ax = nsubplot(1,3,1,1,xsp);
    plot(ax(1),scan10, fastavg(S.Chi(vs,loi),10),'.');
    title(ax(1), sprintf('%s: Line %d', getrun, loi));
    ylabel(ax(1),S.base);
    xlabel(ax(1),'Scan Number');
  else
    ax = [ nsubplot(2,3,1,1,xsp) nsubplot(2,3,2,1,xsp) ];
    plot(ax(1),scan010, fastavg(S0.Chi(vs0,loi),10),'.');
    title(ax(1), sprintf('%s: Line %d', getrun, loi));
    ylabel(ax(1),sprintf('Ref: %s',S0.base));
    set(ax(1),'XTickLabel',[]);
    plot(ax(2),scan10, fastavg(S.Chi(vs,loi),10),'.');
    ylabel(ax(2),sprintf('Ref: %s',S0.base));
    xlabel(ax(2),'Scan Number');
    set(ax(2),'YAxisLocation','Right');
    linkaxes(ax,'x');
  end

  Chimean = mean(S.Chi(vs,loi));
  Chires = (S.Chi(vs,loi)-Chimean)./Chimean;
  % Gedmean = mean(S.Ged(vs,loi));
  Gedrat = S.Ged./S.Gedcalc;
  Gedres = Gedrat(vs,loi);
  fn3 = interp2(XX,YY,fn,line_pos(vs,loi),scans);
  
  Navg = 10;
  fn310 = fastavg(fn3,Navg);
  ffn3 = fn310 - floor(fn310);
  Chires10 = fastavg(Chires,Navg);

  ax = [ nsubplot(2,3,1,2,xsp) nsubplot(2,3,2,2,xsp) ];
  plot(ax(1),scan10,Chires10,'.');
  title(ax(1),sprintf('%s: Line %d', getrun, loi));
  set(ax(1),'XTickLabel',[]);
  ylabel(ax(1),'Mixing Ratio RDFM');
  plot(ax(2),scan10,ffn3,'.');
  set(ax(2),'YAxisLocation','Right');
  xlabel(ax(2),'Scan Number');
  ylabel(ax(2),'Fractional Fringe');
  set(ax(2),'YAxisLocation','Right');
  linkaxes(ax,'x');
  %
  xffn = [0, 0.5, 1];
  vxffn = 0.5*[1, -1, 1];
  vffn = interp1(xffn,vxffn,ffn3);
  Vffn = polyfit(vffn,Chires10,1);
  vffnfit = polyval(Vffn,vxffn);
  osc(loi) = Vffn(1);
  
  ax = nsubplot(2,3,[2 2],3,xsp);
  %scatter(ffn3,Chires10,[],scan10);
  plot(ax,ffn3,Chires10,'.',xffn,vffnfit,'r.-');
  title(ax,sprintf('%s: Line %d', getrun, loi));
  xlabel(ax,'Etalon Fraction fringe position at line center');
  ylabel(ax,'Mixing Ratio RDFM');
end

osc_out = osc;

function [scans1_out, vs] = scan_intersect(scans1, scans2)
if isrow(scans1)
  scans1 = scans1';
end
vi = interp1(scans1,1:length(scans1),scans2,'nearest','extrap');
vs = scans2 == scans1(vi);
scans1_out = scans1(vi(vs));
