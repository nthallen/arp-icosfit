function taus = rawview( dirname, figno );
% taus = rawview( dirname[, figno] );
% Review raw ringdown files.
% Displays ringdown, the onboard linear fit and also a full non-linear
% log fit for comparison.
% Output columns are
%    tau tau2 offset fileno length(raw) std1 std2 k k2
%   tau, std1, k all pertain to the linear fit
%   tau2, std2, k2 pertain to the non-linear least squares
if nargin < 2
  figno = -1;
  if nargin < 1
    dirname = [];
    cr_cfg = load_cr_cfg;
    dirs = { 'RAW', [ cr_cfg.Matlab_CD_Path '/home/CR/' getrun(1) '/CPCI/RAW' ] };
    for i=1:length(dirs)
      if exist( dirs{i}, 'dir' )
        dirname = dirs{i};
        break;
      end
    end
  end
end
% set to 1 if dlmread treats a leading space as an empty column
dlmkluge = 1;
if ~exist( dirname, 'dir' )
  error('Cannot locate RAW directory');
end
% addpath c:/home/analogic
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
      taus = [ taus; rawview( path, figno ) ];
    elseif length(findstr(file.name, '.DAT')) > 0 | length(findstr(file.name, '.dat')) > 0
      try
        rawf = dlmread( [ dirname '/' file.name ], ' ', 0, 0 );
        clf;
        offset = rawf(1,5);
        fileno = rawf(1,6);
        % b = rawf(2,2);
        std1 = sqrt(rawf(1,5+dlmkluge));
        n = rawf(2,6+dlmkluge);
        % tau = n/log(b);
        skip = rawf(2,7+dlmkluge);
        raw = rawf([3:size(rawf,1)],1);
        allv = [1:length(raw)]';
        fitv = [(skip+1):length(raw)]';
        % trialx = exp(-(fitv-skip)/tau);
        % z = a/(1-b);
        % k = sum((raw(fitv)-z).*trialx)./sum(trialx.*trialx);
        % fit = z+k*trialx;
        
        % Do linear auto-correlation fit:
        V = fitlin(raw(fitv), n);
        b = V(2);
        a = V(5);
        tau = n/log(b);
        z = a/(1-b);
        trialx = exp(-(fitv-skip)/tau);
        k = sum((raw(fitv)-z).*trialx)./sum(trialx.*trialx);
        fit = z+k*trialx;
        std1 = sqrt(V(4));
        
        
        % Now do a non-linear logarithmic fit
        V = [ k tau z ];
        % V = fmins('logchi', V, [], [], fitv-skip, raw(fitv) );
        V = fminsearch('logchi', V, [], fitv-skip, raw(fitv) );
        k2 = V(1);
        tau2 = V(2);
        z2 = V(3);
        fit2 = k2*exp(-(fitv-skip)/tau2) + z2;
        std2 = std(raw(fitv)-fit2);
        taus = [ taus; tau tau2 offset fileno length(raw) std1 std2 k k2 ];
        fprintf( 1, '%2d %4d %3d %f %f %f %f %f %f %f %f\n', ...
          offset, fileno, length(raw), ...
          tau, tau2, std1, std2, k, k2, z, z2 );
        
        if figno > -2
          if ishandle(figno)
            figure(figno);
          else
            return
          end
          nsubplot( 3, 1, [ 2 2 ], 1 );
          plot( allv, raw, '.', fitv, fit, fitv, fit2 );
          set(gca,'XTickLabel',[]);
          xl = xlim;
          title( sprintf( '%d: tau1=%f tau2=%f', fileno, tau, V(2) ));
          nsubplot( 3, 1, [ 3 1 ], 1 );
          plot(fitv, raw(fitv)-fit,fitv,raw(fitv)-fit2);
          xlim(xl);
          grid;
          set(gca,'YAxisLocation','right');
          nsubplot(4,4,1.2,3.8);
          hist(raw(fitv)-fit,20);
          drawnow; shg;
          pause;
        end
      catch
        fprintf(1,'Error: %s\n', lasterr);
      end
    end
  end
end
