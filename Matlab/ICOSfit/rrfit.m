function rrfit( base, range, plotcode );
% rrfit( base, [,scans] [, plotcode] );
% Display fitted scan, optionally limited to a range
% of Scan values.
%plotcodes:
%   1: Plots fit and residual
%   2: Plots fit in transmission (baseline removed) and residual
%   3: Plots line strength and residual
%   4,5,6: Same as 1,2,3, but plots versus sample number instead of
%           wavenumber.
%   7,8,9: Same as 1,2,3, but also plots the baseline and detrended
%           baseline on seperate axes.
if nargin < 1
  disp('Error: Must define base. eg ''ICOSout.R1.3p''')
  return
end
ICOSsetup;
if nargin == 1 || isempty(range)
  scans = listscans(base);
else
  if size(range,2) > 1
    range = range';
  end
  rows = unique(interp1( scannum, [1:length(scannum)]', range, 'nearest' ));
  rows = rows(isfinite(rows));
  scans = scannum(rows);
end
if nargin < 3
    plotcode = 1;
end
plotbase = 0;
if plotcode >= 7
    plotbase = 1;
    plotcode = plotcode - 6;
end
if plotcode > 9
    disp('Error: Invalid plotcode.');
    return
end

AppData.base = base;
AppData.plotcode = plotcode;
AppData.plotbase = plotbase;
AppData.scannum = scannum;
AppData.nu = nu;
AppData.nu0 = nu0;
AppData.delta = delta;
AppData.P_vec = P_vec;
AppData.fitdata = fitdata;
AppData.nu_F0 = nu_F0;
AppData.Scorr = Scorr;
AppData.CavLen = CavLen;
AppData.Nfit = Nfit;
AppData.Ged = Ged;
AppData.v = v;

if plotbase == 0
Axes = [
    60    45    60     1    20    15     0     .5
    60    45    60     1     0    45    60     1
    ];
elseif plotbase == 1
    Axes = [
    60    45    60     1    20    15     0     .5
    60    45    60     1     0    45     0     1
    60    45    60     1     0    30     0     .5
    60    45    60     1     0    30    60     .5
    ];
end
scan_viewer('Scans', scans, 'Axes', Axes, 'Name', 'rrfit Viewer', ...
    'Callback', @rrfitview_callback, 'AppData', AppData);

function rrfitview_callback(handles, sv_axes)
if nargin < 2
    sv_axes = handles.Axes;
end
AppData = handles.data.AppData;
scan = handles.data.Scans(handles.data.Index);
path = mlf_path( AppData.base, scan, '.dat');
fe = load( path );
data_ok = (~isempty(fe));
if data_ok
    i=find(AppData.scannum==scan);
    lpos = AppData.nu + AppData.delta*AppData.P_vec(i)/760. - AppData.fitdata(i,AppData.v);
    if AppData.nu0 ~= 0
      lpos = lpos + AppData.nu_F0(i);
    end
    lgd = AppData.fitdata(i,AppData.v+1);
    lst = AppData.Scorr.*AppData.CavLen.*AppData.Nfit./(AppData.Ged*sqrt(pi));
    lgl = AppData.fitdata(i,AppData.v+3);
    lwid = (lgd+lgl); % fitdata(i,v+1)+fitdata(i,v+3);
     if AppData.plotcode <= 3
        nux = fe(:,2);
        if min(nux) < 2
            nux = nux+AppData.nu_F0(i)+AppData.nu0; % This may change
        end
        xdir = 'reverse';
        ttlx = 'Wavenumber (cm-1)';
    elseif AppData.plotcode > 3
        nux = fe(:,1);
        xdir = 'normal';
        ttlx = 'Sample Number';
    end
    dataX = nux;
    resX = nux;
    if AppData.plotcode == 1 || AppData.plotcode == 4
        dataY = fe(:,[3:5]);
        resY = fe(:,3)-fe(:,4);
        reslbl = 'Fit Res';
        sdev = sqrt(mean((resY).^2));
        ttltext = [ 'Scan:' num2str(AppData.scannum(i)) ...
            ' \sigma = ' num2str(sdev)  ' '  path ];
        basep = interp1(nux,fe(:,5),lpos,'nearest');
        fitp =  interp1(nux,fe(:,4),lpos,'nearest');
        meanp = (basep+fitp)/2;
        % threw in maxp to see the position of small lines
        maxp = ones(size(lpos))*max(fe(:,5));
        
    elseif AppData.plotcode == 2 || AppData.plotcode == 5
        dataY = [fe(:,3)./fe(:,5)*100, fe(:,4)./fe(:,5)*100];
        resY = (fe(:,3)-fe(:,4))./fe(:,5)*100;
        reslbl = 'Fit Res (%)';
        sdev = sqrt(mean((resY).^2));
        ttltext = [ 'Scan:' num2str(AppData.scannum(i)) ...
            ' \sigma = ' num2str(sdev)  ' '  path ];
        basep = interp1(nux,fe(:,5)./fe(:,5)*100,lpos,'nearest');
        fitp =  interp1(nux,fe(:,4)./fe(:,5)*100,lpos,'nearest');
        meanp = (basep+fitp)/2;
        % threw in maxp to see the position of small lines
        maxp = ones(size(lpos))*100;
    elseif AppData.plotcode == 3 || AppData.plotcode == 6
        dataY = fe(:,6);
        resY = fe(:,3)-fe(:,4);
        reslbl = 'Fit Res';
        sdev = sqrt(mean((resY).^2));
        ttltext = [ 'Scan:' num2str(AppData.scannum(i)) ...
            ' \sigma = ' num2str(sdev)  ' '  path ];
        basep = interp1(nux,fe(:,6),lpos,'nearest');
        fitp =  interp1(nux,fe(:,6),lpos,'nearest');
        meanp = (basep+fitp)/2;
        % threw in maxp to see the position of small lines
        maxp = ones(size(lpos))*max(fe(:,6));
    end
   
    X = [ lpos lpos-lwid; lpos lpos+lwid ];
    % Y = [ basep meanp; fitp meanp ];
    Y = [ maxp meanp; fitp meanp ];
      
    plot(sv_axes(1),resX, resY);
    set(sv_axes(1),'XDir',xdir,'xticklabel',[],'Xgrid','on','Ygrid','on');
    ylabel(sv_axes(1),reslbl);
    title(sv_axes(1),ttltext);

    if AppData.plotcode <= 3
        plot(sv_axes(2), dataX,dataY, X, Y, 'r');
    elseif AppData.plotcode > 3
        plot(sv_axes(2), dataX,dataY);
    end
    set(sv_axes(2),'XDir', xdir,'YAxisLocation','right','Xgrid','on','Ygrid','on');
    if AppData.plotbase == 1
          set(sv_axes(2),'XTickLabel',[]);
          plot(sv_axes(3),dataX,fe(:,5),'r'); 
          set(sv_axes(3),'XTickLabel',[],'XDir',xdir,'Xgrid','on','Ygrid','on');
          ylabel(sv_axes(3),'Baseline');
          plot(sv_axes(4),dataX,detrend(fe(:,5)),'r');
          set(sv_axes(4),'XDir',xdir,'YAxisLocation','right','Xgrid','on','Ygrid','on');
          ylabel(sv_axes(4),'Detrended Baseline'); 
          xlabel(sv_axes(4),ttlx);
    else
          xlabel(sv_axes(2),ttlx);
    end
end
