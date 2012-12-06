function ringview( scannum, wavenum )
% ringview( [ scannum [, wavenum ]]] );
% Reviews ringdown data.
PT = load('PT');
Waves = load_waves;
cell_cfg=load_cell_cfg;
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
AppData.CavityLength=cell_cfg.CavityLength;
AppData.QCLI_Wave = PT.QCLI_Wave;
AppData.idx = idx;
AppData.wavenum = wavenum;
AppData.Axes = [
    60    45    60     1    20    15    35     .5   0
    60    45    60     1     0    45    50     1    1
    ];

scan_viewer('Scans', scannum, 'Axes', AppData.Axes, 'Name', 'Ringdown Viewer', ...
    'Callback', @ringview_callback, 'AppData', AppData);

function ringview_callback(handles, sv_axes)
if nargin < 2
    sv_axes = handles.Axes;
end
AppData = handles.data.AppData;
if ~isfield(AppData,'menus')
    top_menu = uimenu(handles.figure,'Tag','ringview','Label','ringview');
    cb = @ringview_menu_callback;
    AppData.menus.Fit = uimenu(top_menu,'Tag','fitmenu','Label','Fit Menu');
    AppData.menus.Fit_nonlin = uimenu(AppData.menus.Fit,'Tag','nonlin','Label','Non-linear lsq','Callback',cb);
    AppData.menus.Fit_auto = uimenu(AppData.menus.Fit,'Tag','auto','Label','Autocorrelation','Callback',cb);
    AppData.menus.Fit_both = uimenu(AppData.menus.Fit,'Tag','both','Label','Both','Checked','on','Callback',cb);
    AppData.menus.Tau = uimenu(top_menu,'Tag','taumenu','Label','Tau Display');
    AppData.menus.Tau_current = uimenu(AppData.menus.Tau,'Tag','Tau_current','Label','Current','Callback',cb);
    AppData.menus.Tau_sample = uimenu(AppData.menus.Tau,'Tag','Tau_sample','Label','Sample','Checked','on','Callback',cb);
    AppData.menus.Select = uimenu(top_menu,'Tag','Select','Label','Select Baseline','Callback',cb);
    handles.data.AppData = AppData;
    guidata(handles.figure,handles);
end
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
      v = find(~isnan(fe(:,1)));
      if ~isempty(v)
            xdata=1/AppData.Waves.RawRate*AppData.Waves.NAverage*[1:length(fe(:,1))];
            dt =  mean(diff(xdata));
            n = 5;
            delay = 3.3e-6; %Delay in seconds of the VtoI/electronics
            skip = ceil(AppData.Waves.TzSamples + delay*AppData.Waves.RawRate); %number of points to skip
            xdata=xdata-xdata(skip);  
            if size(fe,1) > skip+1000
                fitv = [skip:1000]';  %points to include in fit
            else
                fitv = [skip:length(fe(:,1))]';
            end
            
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
        text(0.02,0.98, ...
            sprintf('Tau_{auto} = %.2f \\musec (R = %.1f ppm)',nanmedian(AppData.taus(1).Tau)*1e6,AppData.CavityLength/nanmedian(AppData.taus(1).Tau)/2.998e10*1e6), ...
            'Parent',sv_axes(1),'Color','b','Units','Normalized','VerticalAlignment','top');
        text(0.98,0.98, ...
            sprintf('Tau_{nonlin} = %.2f \\musec (R = %.1f ppm)',nanmedian(AppData.taus(2).Tau)*1e6,AppData.CavityLength/nanmedian(AppData.taus(2).Tau)/2.998e10*1e6), ...
            'Parent',sv_axes(1),'Color','g','Units','Normalized','VerticalAlignment','top','HorizontalAlignment','right');
        
        
        plot(sv_axes(2),xdata*1e6,fe(:,1),'k', ...
            [0,0],ylim(sv_axes(2)),':k',xlim(sv_axes(2)),[mean(fe(end-200:end,1)),mean(fe(end-200:end,1))],':k', ...
            xdata(fitv)*1e6,AppData.taus(1).Fit(:,iscan),'b',...
            xdata(fitv)*1e6,AppData.taus(2).Fit(:,iscan),'g');
        xlabel(sv_axes(2),'\musec');
        ylabel(sv_axes(2),'Power');
        text(0.5,0.95, ...
            sprintf('Ringdown Scan: %d\n\nTau_{auto} = %.2f \\musec (std = %.2f)\nTau_{nonlin} = %.2f \\musec (std = %.2f)', ...
            AppData.scannum(iscan), ...
            AppData.taus(1).Tau(iscan)*1e6,AppData.taus(1).Std(iscan), ...
            AppData.taus(2).Tau(iscan)*1e6,AppData.taus(2).Std(iscan) ), ...
            'Parent',sv_axes(2),'Units','Normalized','VerticalAlignment','top','HorizontalAlignment','left');
      end
    end
    
  end

function ringview_menu_callback(hObject,eventdata)
handles = guidata(hObject);
AppData = handles.data.AppData;
Tag = get(hObject,'Tag');
switch Tag(1)
    case 'Y'
        sig = 'off';
        trans = 'off';
        stren = 'off';
        switch hObject
            case AppData.menus.Y_signal
                sig = 'on';
                handles.data.AppData.Yopt = 0;
            case AppData.menus.Y_transmission
                trans = 'on';
                handles.data.AppData.Yopt = 1;
            case AppData.menus.Y_strength
                stren = 'on';
                handles.data.AppData.Yopt = 2;
        end
        set(AppData.menus.Y_signal,'checked',sig);
        set(AppData.menus.Y_transmission,'checked',trans);
        set(AppData.menus.Y_strength,'checked',stren);
        handles.data.ylim{2} = [];
        guidata(hObject,handles);
        scan_viewer('scan_display',handles);
    case 'X'
        wvno = 'off';
        samp = 'off';
        switch hObject
            case AppData.menus.X_wavenumber
                wvno = 'on';
                handles.data.AppData.Xopt = 0;
            case AppData.menus.X_sample
                samp = 'on';
                handles.data.AppData.Xopt = 1;
        end
        for i = 1:length(handles.data.xlim)
            handles.data.xlim{i} = [];
        end
        set(AppData.menus.X_wavenumber,'checked',wvno);
        set(AppData.menus.X_sample,'checked',samp);
        guidata(hObject,handles);
        scan_viewer('scan_display',handles);
    case 'B'
        if strcmp(get(AppData.menus.Baselines,'checked'),'on')
            set(AppData.menus.Baselines,'checked','off');
            handles.data.AppData.plotbase = 0;
            handles.data.Axes = AppData.Axes_A;
        else
            set(AppData.menus.Baselines,'checked','on');
            handles.data.AppData.plotbase = 1;
            handles.data.Axes = AppData.Axes_B;
        end
        guidata(hObject,handles);
        scan_viewer('figure_ResizeFcn',hObject,eventdata,handles);
        handles = guidata(hObject);
        scan_viewer('scan_display',handles);
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
