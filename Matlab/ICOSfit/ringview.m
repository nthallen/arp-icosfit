function ringview( scannum, wavenum, basepath)
% ringview( [ scannum [, wavenum[, basepath]]]] );
% Reviews ringdown data.
if nargin < 3
    basepath = '';
end
PT = load('PT');
[Waves,WaveRange] = waves_used;
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
wavenumi = struct2cell(WaveRange);
wavenumi = cell2mat(squeeze(wavenumi(1,:,:)));

roris = ~[ Waves(wavenums==wavenumi).ISICOS ];
if nargin >= 2 && ~isempty(wavenum)
  roris = roris & (wavenums == wavenum)';
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

AppData.WaveRange = WaveRange(wavenum==wavenumi).ranges;
iring=find(PT.QCLI_Wave(idx)==wavenum);
scannum=scannum(iring);
idx=idx(iring);
AppData.Waves = Waves(wavenum==wavenumi);
AppData.base = find_scans_dir(basepath);
AppData.binary = 1;
if size(scannum,1) > 1; scannum = scannum'; end
AppData.scannum = scannum;
AppData.CavityLength=cell_cfg.CavityLength;
AppData.QCLI_Wave = PT.QCLI_Wave;
path = mlf_path( AppData.base, AppData.WaveRange(1), '.dat');
[fe, hdr] = loadbin(path);
AppData.StartSerNum = hdr.SerialNum;
AppData.idx = idx;
AppData.wavenum = wavenum;
%Setup x vector in microsec, correlation shift, and delay.
xdata=1/AppData.Waves.RawRate*AppData.Waves.NAverage*[1:length(fe(:,1))];
AppData.dt =  mean(diff(xdata));
AppData.n = 2; %Correlation shift
AppData.delay = 3.5e-6; %Delay in seconds of the VtoI/electronics
AppData.skip = ceil(AppData.Waves.TzSamples + AppData.delay*AppData.Waves.RawRate); %number of points to skip
AppData.xdata=xdata-xdata(AppData.skip);  
if size(fe,1) > AppData.skip+1000
    AppData.fitv = [AppData.skip:1000]';  %points to include in fit
else
    AppData.fitv = [AppData.skip:length(fe(:,1))]';
end
AppData.FitDisplay = 2;
AppData.TauDisplay = 1;
AppData.tauwindow.x = [];
AppData.tauwindow.y = [];
AppData.taus = struct('Name',{'auto','nonlin'}, ...
    'Tau',{zeros(size(scannum))*NaN,zeros(size(scannum))*NaN}, ...
    'Std',{zeros(size(scannum))*NaN,zeros(size(scannum))*NaN}, ...
    'Fit',{zeros(length(AppData.fitv),length(scannum))*NaN,zeros(length(AppData.fitv),length(scannum))*NaN}, ...
    'MeanTau',{[],[]}, ...
    'ScanNum',{scannum,scannum}, ...
    'SerialNum',{zeros(size(scannum))*NaN}, ...
    'CurrentOffset',{zeros(size(scannum))*NaN}, ...
    'Etalon',{zeros(size(scannum))*NaN}, ...
    'Status',{zeros(size(scannum))*NaN} );
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
    top_menu = uimenu(handles.scan_viewer,'Tag','ringview','Label','ringview');
    cb = @ringview_menu_callback;
    AppData.menus.Fit = uimenu(top_menu,'Tag','fitmenu','Label','Fit Menu');
    AppData.menus.Fit_nonlin = uimenu(AppData.menus.Fit,'Tag','Fit_nonlin','Label','Non-linear lsq','Callback',cb);
    AppData.menus.Fit_auto = uimenu(AppData.menus.Fit,'Tag','Fit_auto','Label','Autocorrelation','Callback',cb);
    AppData.menus.Fit_both = uimenu(AppData.menus.Fit,'Tag','Fit_both','Label','Both','Checked','on','Callback',cb);
    AppData.menus.Tau = uimenu(top_menu,'Tag','taumenu','Label','Tau Display');
    if AppData.Waves.NetSamples > 1
    AppData.menus.Tau_current = uimenu(AppData.menus.Tau,'Tag','Tau_current','Label','Current','Callback',cb);
    end
    AppData.menus.Tau_sample = uimenu(AppData.menus.Tau,'Tag','Tau_sample','Label','Sample','Checked','on','Callback',cb);
    AppData.menus.Export = uimenu(top_menu,'Tag','Export','Label','Export Tau Struct','Callback',cb);
    AppData.menus.Select = uimenu(top_menu,'Tag','Select','Label','Select Base Tau Region','Callback',cb);
    AppData.menus.Write = uimenu(top_menu,'Tag','Write','Label','Write to Cell_Config','Callback',cb);
    handles.data.AppData = AppData;
    guidata(handles.scan_viewer,handles);
end
scan = handles.data.Scans(handles.data.Index); %scan number
iscan = find(AppData.scannum == scan); %index for scan into scannum and idx
if AppData.QCLI_Wave(AppData.idx(iscan)) == AppData.wavenum
    path = mlf_path( AppData.base, scan, '.dat');
    data_ok = 0;
    if AppData.binary
      [fe,hdr] = loadbin( path );
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
          
        %Fit data if not already fit    
        if isnan(AppData.taus(1).Tau(iscan)) && hdr.SerialNum ~= AppData.StartSerNum
            AppData.taus(1).CurrentOffset(iscan) = mod(hdr.SerialNum-AppData.StartSerNum,AppData.Waves.NetSamples);
            AppData.taus(2).CurrentOffset(iscan) = mod(hdr.SerialNum-AppData.StartSerNum,AppData.Waves.NetSamples);

            AppData.taus(1).SerialNum(iscan) = hdr.SerialNum;
            AppData.taus(2).SerialNum(iscan) = hdr.SerialNum;

            AppData.taus(1).Etalon(iscan) = fe(50,2) - min(fe(:,2));
            AppData.taus(2).Etalon(iscan) = fe(50,2) - min(fe(:,2));

            AppData.taus(1).Status = hdr.Status;
            AppData.taus(2).Status = hdr.Status;

            %Do linear auto-correlation fit:
            V = fitlin(fe(AppData.fitv,1), AppData.n);
            b = V(2);
            a = V(5);
            tau = AppData.n*AppData.dt/log(b);
            z = a/(1-b);
            trialx = exp(-AppData.xdata(AppData.fitv)/tau);
            k = sum((fe(AppData.fitv,1)'-z).*trialx)./sum(trialx.*trialx);
            fit = z+k*trialx;
            std1 = sqrt(V(4));
            AppData.taus(1).Tau(iscan) = tau;
            AppData.taus(1).Std(iscan) = std1;
            AppData.taus(1).Fit(:,iscan) = fit;
            
            % Now do a non-linear logarithmic fit
            V = [ k tau z ];
            V = fminsearch('logchi', V, [], AppData.xdata(AppData.fitv)', fe(AppData.fitv,1) );
            k2 = V(1);
            tau2 = V(2);
            z2 = V(3);
            fit2 = k2*exp(-AppData.xdata(AppData.fitv)/tau2) + z2;
            std2 = std(fe(AppData.fitv,1)'-fit2);
            AppData.taus(2).Tau(iscan) = tau2;
            AppData.taus(2).Std(iscan) = std2;
            AppData.taus(2).Fit(:,iscan) = fit2;
            handles.data.AppData = AppData;
            guidata(handles.scan_viewer,handles)
        end
        %Set x units for tau display 
        if AppData.TauDisplay == 1
            xlab = 'Scan Number';
            xtau = AppData.taus.ScanNum;
        elseif AppData.TauDisplay == 0
            xlab = 'Current Number';
            xtau = AppData.taus.CurrentOffset;
        end     
        %Plot fitted taus
        newplot(sv_axes(1))
        xlabel(sv_axes(1),xlab)
        ylabel(sv_axes(1),'Tau (\musec)')
        title(sv_axes(1),getrun(0,handles.scan_viewer))
        if AppData.FitDisplay == 1 || AppData.FitDisplay == 2
            line(xtau,AppData.taus(1).Tau*1e6,'Parent',sv_axes(1),'Color','b','LineStyle','none','Marker','.')
            text(0.02,0.98, ...
                sprintf('Tau_{auto} = %.2f \\musec (R = %.1f ppm)',nanmedian(AppData.taus(1).Tau)*1e6,AppData.CavityLength/nanmedian(AppData.taus(1).Tau)/2.998e10*1e6), ...
                'Parent',sv_axes(1),'Color','b','Units','Normalized','VerticalAlignment','top');
        end
        if AppData.FitDisplay == 0 || AppData.FitDisplay == 2
            line(xtau,AppData.taus(2).Tau*1e6,'Parent',sv_axes(1),'Color','g','LineStyle','none','Marker','.')
            text(0.98,0.98, ...
                sprintf('Tau_{nonlin} = %.2f \\musec (R = %.1f ppm)',nanmedian(AppData.taus(2).Tau)*1e6,AppData.CavityLength/nanmedian(AppData.taus(2).Tau)/2.998e10*1e6), ...
                'Parent',sv_axes(1),'Color','g','Units','Normalized','VerticalAlignment','top','HorizontalAlignment','right');
        end
        if isempty(handles.data.xlim{1})
            xlim(sv_axes(1),'auto')
        end
        yl=ylim(sv_axes(1));
        if yl(1) < 0; yl(1) = 0; end
        if yl(2) <= yl(1); yl(2) = yl(1) + 1; end
        if yl(2) > 50; yl(2) = 50; end
        ylim(sv_axes(1),yl);
        %Calculate mean tau in select box region
        if ~isempty(AppData.tauwindow.x)
           line(AppData.tauwindow.x,AppData.tauwindow.y,'Color','k','LineWidth',2,'Parent',sv_axes(1)) 
           x = AppData.tauwindow.x;
           y = AppData.tauwindow.y;
           itauauto = xtau>=min(x) & xtau<=max(x) & AppData.taus(1).Tau>=min(y)*1e-6 & AppData.taus(1).Tau<=max(y)*1e-6;
           itaunonlin = xtau>=min(x) & xtau<=max(x) & AppData.taus(2).Tau>=min(y)*1e-6 & AppData.taus(2).Tau<=max(y)*1e-6;
           AppData.taus(1).MeanTau = nanmean(AppData.taus(1).Tau(itauauto));
           AppData.taus(2).MeanTau = nanmean(AppData.taus(2).Tau(itaunonlin));
           handles.data.AppData = AppData;
           guidata(handles.scan_viewer,handles);
        end
        %Display mean tau
        if ~isempty(AppData.taus(1).MeanTau)
            if AppData.FitDisplay == 1 || AppData.FitDisplay == 2
                text(0.02,0.85, ...
                    sprintf('Tau_{auto} = %.2f \\musec (R = %.1f ppm)',AppData.taus(1).MeanTau*1e6,AppData.CavityLength/AppData.taus(1).MeanTau/2.998e10*1e6), ...
                    'Parent',sv_axes(1),'Units','Normalized','VerticalAlignment','top');
            end
            if AppData.FitDisplay == 0 || AppData.FitDisplay == 2
                text(0.98,0.85, ...
                    sprintf('Tau_{nonlin} = %.2f \\musec (R = %.1f ppm)',AppData.taus(2).MeanTau*1e6,AppData.CavityLength/AppData.taus(2).MeanTau/2.998e10*1e6), ...
                    'Parent',sv_axes(1),'Units','Normalized','VerticalAlignment','top','HorizontalAlignment','right');
            end 
        end
        %Plot individual fits
        newplot(sv_axes(2));
        plot(sv_axes(2),AppData.xdata*1e6,fe(:,1),'k');
        line([0,0],ylim(sv_axes(2)),'Color','k','LineStyle',':','Parent',sv_axes(2));
        line(xlim(sv_axes(2)),[mean(fe(end-200:end,1)),mean(fe(end-200:end,1))],'Color','k','LineStyle',':','Parent',sv_axes(2));
        if AppData.FitDisplay == 1 || AppData.FitDisplay == 2
            line(AppData.xdata(AppData.fitv)*1e6,AppData.taus(1).Fit(:,iscan),'Parent',sv_axes(2),'Color','b')
        end
        if AppData.FitDisplay == 0 || AppData.FitDisplay == 2
            line(AppData.xdata(AppData.fitv)*1e6,AppData.taus(2).Fit(:,iscan),'Parent',sv_axes(2),'Color','g')
        end
        xlabel(sv_axes(2),'\musec');
        ylabel(sv_axes(2),'Power');
        tautext = sprintf('Ringdown Scan: %d', AppData.scannum(iscan));
        text(0.5,0.95, tautext, ...
            'Parent',sv_axes(2),'Units','Normalized','VerticalAlignment','top','HorizontalAlignment','left');      
        if AppData.FitDisplay == 0 || AppData.FitDisplay == 2
            tautext = sprintf('Tau_{nonlin} = %.2f \\musec (std = %.2f)', ...
            AppData.taus(2).Tau(iscan)*1e6,AppData.taus(2).Std(iscan) );
            text(0.5,0.8, tautext, ...
                'Parent',sv_axes(2),'Color','g','Units','Normalized','VerticalAlignment','top','HorizontalAlignment','left');      
        end
        if AppData.FitDisplay == 1 || AppData.FitDisplay == 2
            tautext = sprintf('Tau_{auto} = %.2f \\musec (std = %.2f)', ...
            AppData.taus(1).Tau(iscan)*1e6,AppData.taus(1).Std(iscan) );
            text(0.5,0.7, tautext, ...
                'Parent',sv_axes(2),'Color','b','Units','Normalized','VerticalAlignment','top','HorizontalAlignment','left');
        end
     end
   end  
end

function ringview_menu_callback(hObject,eventdata)
handles = guidata(hObject);
AppData = handles.data.AppData;
Tag = get(hObject,'Tag');
switch Tag(1)
    case 'F'
        nonlin = 'off';
        auto = 'off';
        both = 'off';
        switch hObject
            case AppData.menus.Fit_nonlin
                nonlin = 'on';
                AppData.FitDisplay = 0;
            case AppData.menus.Fit_auto
                auto = 'on';
                AppData.FitDisplay = 1;
            case AppData.menus.Fit_both
                both = 'on';
                AppData.FitDisplay = 2;
        end
        set(AppData.menus.Fit_nonlin,'Checked',nonlin);
        set(AppData.menus.Fit_auto,'Checked',auto);
        set(AppData.menus.Fit_both,'Checked',both);
        handles.data.AppData = AppData;
        guidata(hObject,handles);
        scan_viewer('scan_display',handles);
    case 'T'
        current = 'off';
        sample = 'off';
        AppData.tauwindow.x = [];
        AppData.tauwindow.y = [];
        handles.data.xlim{1} = [];
        switch hObject
            case AppData.menus.Tau_current
                current = 'on';
                AppData.TauDisplay = 0;
            case AppData.menus.Tau_sample
                sample = 'on';
                AppData.TauDisplay = 1;
        end
        set(AppData.menus.Tau_current,'checked',current);
        set(AppData.menus.Tau_sample,'checked',sample);
        handles.data.AppData = AppData;
        guidata(hObject,handles);
        scan_viewer('scan_display',handles);
    case 'S'
        k = waitforbuttonpress;
        point1 = get(gca,'CurrentPoint');    % button down detected
        finalRect = rbbox;                   % return figure units
        point2 = get(gca,'CurrentPoint');    % button up detected
        point1 = point1(1,1:2);              % extract x and y
        point2 = point2(1,1:2);
        p1 = min(point1,point2);             % calculate locations
        offset = abs(point1-point2);         % and dimensions
        AppData.tauwindow.x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
        AppData.tauwindow.y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
        handles.data.AppData = AppData;
        guidata(hObject,handles);
        handles = guidata(hObject);
        scan_viewer('scan_display',handles);
    case 'E'
        assignin('base','tau',AppData.taus);
    case 'W'
        cell_cfg = load_cell_cfg;
        fd = fopen([ 'Cell_Config.m'], 'w');
        fprintf(fd, 'function cell_cfg = Cell_Config;\n');
        fprintf(fd, 'cell_cfg.fsr = %.6f;\n', cell_cfg.fsr );
        fprintf(fd, 'cell_cfg.CavityLength = %.2f;\n', cell_cfg.CavityLength );
        if AppData.FitDisplay == 0
            fprintf(fd, 'cell_cfg.MirrorLoss = %.2f;\n', AppData.CavityLength/AppData.taus(2).MeanTau/2.998e10*1e6 );
        elseif AppData.FitDisplay == 1 || AppData.FitDisplay == 2
           fprintf(fd, 'cell_cfg.MirrorLoss = %.2f;\n', AppData.CavityLength/AppData.taus(1).MeanTau/2.998e10*1e6 );
        end 
        fprintf(fd, 'cell_cfg.N_Passes = %i;\n', cell_cfg.N_Passes );
        fprintf(fd, 'cell_cfg.CavityFixedLength = %.2f;\n', cell_cfg.CavityFixedLength );
        fclose(fd);
end
