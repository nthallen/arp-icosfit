function taus = rawview2( dirname, figno );
% rawview2( dirname[, figno] );
% Review raw ringdown files.
% Attempts a non-linear fit to two exponentials.
if nargin < 2
  figno = -1;
  if nargin < 1
    dirname = [];
    dirs = { 'RAW', ['cd/home/CR/' getrun(1) '/CPCI/RAW' ] };
    for i=1:length(dirs)
      if exist( dirs{i}, 'dir' )
        dirname = dirs{i};
        break;
      end
    end
  end
end
% set to 1 if dlmread treats a leading space as an empty column
dlmkluge = 0;
revlevel = sscanf(version,'%*s (R%d)');
if  revlevel < 13
  dlmkluge = 1;
end
if ~exist( dirname, 'dir' )
  error('Cannot locate RAW directory');
end
taus = [];
files = dir( dirname );
[ sfiles ifiles ] = sort({files.name});
files = files(ifiles);
lasterr('');
if figno == -1
  figno = figure;
end
for file = files'
  if figno > 0 & ~ishandle(figno)
    return
  end
  if file.name(1) ~= '.'
    path = [ dirname '/' file.name ];
    if file.isdir
      taus = [ taus; rawview2( path, figno ) ];
    elseif length(findstr(file.name, '.DAT')) > 0 | length(findstr(file.name, '.dat')) > 0
      rawf = dlmread( [ dirname '/' file.name ], ' ', 0, 0 );
      clf;
      offset = rawf(1,5);
      fileno = rawf(1,6);
      b = rawf(2,2);
      std1 = sqrt(rawf(1+dlmkluge,5));
      n = rawf(2,6+dlmkluge);
      tau = 1e-7*n/log(b);
      skip = rawf(2,7+dlmkluge);
      raw = rawf([3:size(rawf,1)],1);
      raw = (raw - mean(raw))/std(raw);
      allv = [1:length(raw)]';
      fitv = [(skip+1):length(raw)]';

      % Do a non-linear dual-exponential fit
      % Initial guesses:
      t1 = 14e-6;
      t2 = 12.25e-6;
      b = min(raw);
      a = (max(raw)-b);
      V = [ a b t1 t2 ];
      t = [0:length(raw)-1]'*1e-7;
      if revlevel < 13
        V = fmins('twoexp', V, [], [], t, raw );
      else
        V = fminsearch(@twoexp, V, ...
          optimset('MaxIter',10000,'MaxFunEvals',10000), t, raw);
      end
      fit1 = calctwoexp(V,t,raw);
      std1 = std(raw-fit1);
      taus = [ taus; fileno V std1 ];
      b = V(2);
      
      if figno > -2
       if ishandle(figno)
         figure(figno);
       else
         return
       end
       nsubplot( 3, 1, [ 2 2 ], 1 );
       plot( allv, raw-b, '.', allv, fit1-b );
       set(gca,'XTickLabel',[]);
       xl = xlim;
       title( sprintf( '%d: tau1=%.1f usecs tau2=%.1f usecs', ...
         fileno, V(3)*1e6, V(4)*1e6 ));
       nsubplot( 3, 1, [ 3 1 ], 1 );
       plot(allv, raw-fit1, '.');
       xlim(xl);
       grid;
       set(gca,'YAxisLocation','right');
       nsubplot(4,4,1.2,3.8);
       hist(raw-fit1,20);
       drawnow; shg;
	   pause;
      end
    end
  end
end

function fit = calctwoexp( V, t, raw );
a = V(1); b = V(2); t1 = V(3); t2 = V(4);
% fit = a * ( t1*exp(-t/t1) - t2*exp(-t/t2)) + b;
fit = a * exp(-t/t1).*(1 + t/t1 - t.*t*(t1-t2)/(2*t2*t1*t1)) + b;
return;

function chi2 = twoexp( V, t, raw );
global neval chi2s
fit = calctwoexp(V,t,raw);
chi2 = sum((fit-raw).^2);
if 0
  clf;
  plot(t,fit-raw);
  title(sprintf('chi2=%g', chi2));
  pause;
end
return;
