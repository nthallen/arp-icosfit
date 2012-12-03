function ringview( scannum, wavenum )
% ringview( [ scannum [, wavenum ]]] );
% Reviews ringdown data.
PT = load('PT');
Waves = load_waves;
v = find(diff(PT.ScanNum))+1; % index of new ScanNum values
if nargin < 1
  scannum = [];
end
if isempty(scannum)
  scannum = PT.ScanNum;
  scannum = [min(scannum(scannum > 0)):max(scannum)];
else
  % constrain to within the range of defined values
  scannum = scannum(scannum > 0 & scannum >= min(PT.ScanNum) & ...
    scannum <= max(PT.ScanNum));
end
base = '';
% Now locate each specified scannum as an index into v, the
% unique PT.ScanNum entry indexes
idx = v(ceil(interp1( PT.ScanNum(v), [1:length(v)], scannum )));
% idx is now an array as long as scannum
% PT.ScanNum(idx) should be equal to scannum except where skipping
% occurs, and then it should be greater than scannum.
wavenums = unique(PT.QCLI_Wave(idx));
roris = ~[ Waves(wavenums+1).ISICOS ];
if nargin >= 2
  roris = roris & (wavenums == wavenum);
else
% now find(roris) has the ringdown entries
% scannum(find(roris)) are the ones we want
% PT.QCLI_Wave(idx(find(roris))) is the wavenum
  wavenum = wavenums(find(roris));
  if length(wavenum) > 1
    wavenum
    error('More than one ringdown waveform implicated: choose one');
  elseif isempty(wavenum)
    error('No ringdown waveforms found');
  end
end
iring=find(PT.QCLI_Wave(idx)==wavenum);
scannum=scannum(iring);
taus = struct('Name',{'auto','nonlin'}, ...
    'Tau',{zeros(size(scannum))*NaN,zeros(size(scannum))*NaN}, ...
    'Std',{zeros(size(scannum))*NaN,zeros(size(scannum))*NaN}, ...
    'Fit',{[],[]});
AppData.taus = taus;
idx=idx(iring);
AppData.Waves = Waves(wavenum+1);
AppData.base = find_scans_dir('');
AppData.binary = 1;
if size(scannum,1) > 1; scannum = scannum'; end
AppData.scannum = scannum;
AppData.n_currents = Waves(wavenum+1).RawSamples;
AppData.RawRate = Waves(wavenum+1).RawRate;
AppData.common_n = 0;
% if n_currents > 1
%   figno = figure;
% else
%   figno = 0;
% end
  
AppData.QCLI_Wave = PT.QCLI_Wave;
AppData.idx = idx;
AppData.wavenum = wavenum;
AppData.Axes = [
    60    45    60     1    20    15    35     .5   0
    60    45    60     1     0    45    50     1    0
    ];

scan_viewer('Scans', scannum, 'Axes', AppData.Axes, 'Name', 'Ringdown Viewer', ...
    'Callback', @ringview_callback, 'AppData', AppData);

function ringview_callback(handles, sv_axes)
if nargin < 2
    sv_axes = handles.Axes;
end
AppData = handles.data.AppData;
scan = handles.data.Scans(handles.data.Index); %scan number
iscan = find(AppData.scannum == scan); %index for scan into scannum and idx
if AppData.QCLI_Wave(AppData.idx(iscan)) == AppData.wavenum
    path = mlf_path( AppData.base, scan, '.dat');
    data_ok = 0;
    if AppData.binary
      fe = loadbin( path );
      data_ok = 1;
    else
      if exist(path,'file')
        fe = load(path);
        data_ok = 1;
      end
    end
    if data_ok
      nsamples = size(fe,1);
      if nsamples ~= AppData.n_currents
        error('nsample ~= n_currents');
      end
      v = find(~isnan(fe(:,1)));
      if ~isempty(v)
        if AppData.common_n == 0
          AppData.common_n = fe(v(1),1);
        end
          xdata=1/AppData.Waves.RawRate*AppData.Waves.NAverage*[1:length(fe(:,1))];
          xdata=xdata-xdata(100);
            dt =  mean(diff(xdata));
            n = 1;
            skip = 105; %number of points to skip
            %fitv = [(skip+1):length(fe(:,1))]';
            fitv = [(skip+1):1000]';  %points to include in fit
        if isnan(AppData.taus(1).Tau(iscan))
            %Do linear auto-correlation fit:
            V = fitlin(fe(fitv,1), n);
            b = V(2);
            a = V(5);
            tau = n*dt/log(b);
            z = a/(1-b);
            trialx = exp(-xdata(fitv)/tau);
            k = sum((fe(fitv,1)'-z).*trialx)./sum(trialx.*trialx);
            fit = z+k*trialx;
            std1 = sqrt(V(4));
            AppData.taus(1).Tau(iscan) = tau;
            AppData.taus(1).Std(iscan) = std1;
            AppData.taus(1).Fit(:,iscan) = fit;
            
            % Now do a non-linear logarithmic fit
            V = [ k tau z ];
            V = fminsearch('logchi', V, [], xdata(fitv)', fe(fitv,1) );
            k2 = V(1);
            tau2 = V(2);
            z2 = V(3);
            fit2 = k2*exp(-xdata(fitv)/tau2) + z2;
            std2 = std(fe(fitv,1)'-fit2);
            AppData.taus(2).Tau(iscan) = tau2;
            AppData.taus(2).Std(iscan) = std2;
            AppData.taus(2).Fit(:,iscan) = fit2;
            handles.data.AppData = AppData;
            guidata(handles.figure,handles)
        end
        plot(sv_axes(1),AppData.scannum,AppData.taus(1).Tau*1e6,'.b',AppData.scannum,AppData.taus(2).Tau*1e6,'.g')
        %set(sv_axes(1),'YAxisLocation','right')
        xlabel(sv_axes(1),'Scan Number')
        ylabel(sv_axes(1),'Tau (\musec)')
        title(sv_axes(1),getrun)
        text(mean(xlim(sv_axes(1))),mean([max(ylim(sv_axes(1))),nanmean(AppData.taus(1).Tau*1e6)]), ...
            sprintf('Tau_{auto} = %.2f \\musec',nanmean(AppData.taus(1).Tau)*1e6), ...
            'Parent',sv_axes(1),'Color','b');
        text(mean(xlim(sv_axes(1))),mean([min(ylim(sv_axes(1))),nanmean(AppData.taus(2).Tau*1e6)]), ...
            sprintf('Tau_{auto} = %.2f \\musec',nanmean(AppData.taus(2).Tau)*1e6), ...
            'Parent',sv_axes(1),'Color','g');
        
        
        plot(sv_axes(2),xdata*1e6,fe(:,1),'k', ...
            [0,0],ylim(sv_axes(2)),':k',xlim(sv_axes(2)),[mean(fe(end-200:end,1)),mean(fe(end-200:end,1))],':k', ...
            xdata(fitv)*1e6,AppData.taus(1).Fit(:,iscan),'b',...
            xdata(fitv)*1e6,AppData.taus(2).Fit(:,iscan),'g');
        xlabel(sv_axes(2),'\musec');
        ylabel(sv_axes(2),'Power');
        text(mean(xlim(sv_axes(2))),max(ylim(sv_axes(2)))-4e3, ...
            sprintf('Ringdown Scan: %d\n\nTau_{auto} = %.2f \\musec (std = %.2f)\nTau_{nonlin} = %.2f \\musec (std = %.2f)', ...
            AppData.scannum(iscan), ...
            AppData.taus(1).Tau(iscan)*1e6,AppData.taus(1).Std(iscan), ...
            AppData.taus(2).Tau(iscan)*1e6,AppData.taus(2).Std(iscan) ), ...
            'Parent',sv_axes(2));
      end
    end
    
  end


% if any(nbinned > 0)
%   ringbins = ringbins./nbinned;
%   taubin = (1e6/RawRate) * common_n ./ log(ringbins);
%   if n_currents > 1
%     if figno > 0
%       figure(figno);
%     else
%       figure;
%     end
%     x = [1:length(taubin)];
%     plot( x, taubin );
%     title(sprintf('Binned Ringdown Tau: %s', getrun));
%     xlabel('Offset');
%     ylabel('Tau \mu secs');
%     cell_cfg=load_cell_cfg;
%     CavityLength = cell_cfg.CavityLength;
%     fprintf(1,'Select Range to average\n');
%     k = waitforbuttonpress;
%     point1 = get(gca,'CurrentPoint');
%     finalRect = rbbox;
%     point2 = get(gca,'CurrentPoint');
%     point1 = point1(1,1:2);              % extract x and y
%     point2 = point2(1,1:2);
%     p1 = min(point1,point2);             % calculate locations
%     offset = abs(point1-point2);         % and dimensions
%     v = find(x >= p1(1) & x <= p1(1)+offset(1));
%     if length(v) > 0
%       taumean = mean(taubin(v));
%       hold on;
%       plot( [p1(1) p1(1)+offset(1)], [taumean taumean], 'r');
%       hold off;
%       fprintf(1, 'Mean tau is %.2f usecs\n', taumean );
%       c = 2.99792458e10; % cm/s
%       MirrorLoss = CavityLength/(c * taumean * 1e-6);
%       fprintf(1, 'MirrorLoss is %.1f ppm\n', MirrorLoss * 1e6 );
%       save MirrorLoss.mat MirrorLoss
%     end
%   end
% end
% 
% if nargout >= 1
%   Taur = taubin;
% end
