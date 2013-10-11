function [maxabsdev, stddev, maxval] = test_icosfit_recalc(base, varargin)
% [maxabsdev, stddev] = test_icosfit_recalc(base)
% test_icosfit_recalc(base, 'view')

% Tests:
%   Background: Verify that BackgroundRegion is 5:TzSamples
%   Raw: Verify that raw scan and raw column of verbose output
%     are identical
%   Base: Verify that recalculated baseline matches verbose baseline
%   nu: Verify that nu gets calculated correctly
%   Fit: Verify that recalculated fit matches verbose fit
AppData.RC = icosfit_recalculate(base);
S = AppData.RC.S;
if S.ICOSfit_format_ver < 2
    error('ICOSfit recalculation requires ICOSfit_format_ver >= 2');
end
S.nu_P = (S.col*S.nu - S.nu0) + (S.col*S.delta).*S.P/760.;
AppData.ICOSfit_cfg = load_ICOSfit_cfg;
AppData.ScanBase = find_scans_dir([]);
AppData.Scans = S.scannum;
AppData.IScans = 1:length(AppData.Scans);
% now AppData.Scans == S.scannum(AppData.IScans)
% When scans are specified explicitly, this allows us to deal with
% fits that might run in reverse or even have multiple overlapping
% regions.
AppData.S = S;

% Load PTE stuff
if ~exist(AppData.S.PTEfile,'file')
    error('Unable to locate PTEfile %s', AppData.S.PTEfile);
end
AppData.PTE = load(AppData.S.PTEfile);

AppData.opt_view = 0; % Show results interactively
AppData.Xopt = 0;
AppData.Yopt = 4;
AppData.Detrend = 0;
AppData.ylbl = 'Fit';
AppData.plotbase = 0;
i = 1;
while i <= length(varargin)
    if ~ischar(varargin{i})
        error('Expected option string at argument %d', i+1);
    elseif strcmpi(varargin{i}, 'view')
        AppData.opt_view = 1;
    elseif strcmpi(varargin{i}, 'scans')
        i = i + 1;
        if i > length(varargin)
            error('Scans option requires value');
        end
        AppData.Scans = varargin{i};
        if size(AppData.Scans,1) == 1
            AppData.Scans = AppData.Scans';
        end
        % should check that these scans are all in S.scannum
        [ss,I] = unique(AppData.S.scannum);
        Is = interp1(ss,I,AppData.Scans,'nearest','extrap');
        v = AppData.S.scannum(Is) == AppData.Scans;
        if ~any(v)
            error('Specified scan range does not match fit data range');
        elseif ~all(v)
            warning('ICOSfit:recalc:ScanNumRange', 'Some scans in specified range are missing');
            AppData.Scans = AppData.Scans(v);
            AppData.IScans = Is(v);
        else
            AppData.IScans = Is;
        end
    else
        error('Unrecognized option: "%s"', varargin{i});
    end
    i = i + 1;
end

% Now map AppData.Scans onto rows of PTE
[ps,I] = unique(AppData.PTE(:,1));
Ip = interp1(ps,I,AppData.Scans,'nearest','extrap');
v = AppData.PTE(Ip,1) == AppData.Scans;
if ~all(v)
    error('Specified scan range does not match PTE file range');
else
    AppData.PTEi = Ip;
end

wvs = waves_used(AppData.Scans);
if length(wvs) > 1
    error('More than one waveform for specified scans');
end
if any(AppData.S.BackgroundRegion ~= [5 wvs.TzSamples])
    warning('ICOSfit:recalc:BackgroundRegion', 'BackgroundRegion differs from default');
end
AppData.BGi = AppData.S.BackgroundRegion(1):AppData.S.BackgroundRegion(2);
if AppData.opt_view
    scan_viewer('Scans', AppData.Scans, ...
        'Axes', [
        60    45    60     1    20    15     0     .5   0
        60    45    60     1     0    45    60     1    0 ], ...
        'Name', 'Test ICOSfit Recalc Viewer', ...
        'Callback', @recalcview_callback, 'AppData', AppData);
else
    [maxabsdev_i, stddev_i, maxval_i] = blind_test(AppData);
    if nargout >= 3
        maxabsdev = maxabsdev_i;
        stddev = stddev_i;
        maxval = maxval_i;
    end
    figure; nsubplot(2,1,1);
    plot(AppData.Scans, maxabsdev_i);
    set(gca,'xticklabel',[]);
    ylabel('max abs dev');
    nsubplot(2,1,2);
    plot(AppData.Scans, stddev_i);
    set(gca,'yaxislocation','right');
    ylabel('stddev');
    legend('sample','nu','raw','fit','baseline','abs');
end

% function Baseline = loadetlnbase(AppData)
% % Assign Baseline params to structure
% [ nu, vectors, Pdegree, Ptype, PV, Pscale ] = ...
%     readetlnbase( AppData.S.BaselineFile );
% Baseline = struct('nu', nu, 'vectors', vectors, ...
%     'Pdegree', Pdegree, 'Ptype', Ptype, 'PV', PV, ...
%     'Pscale', Pscale);
% assert(Ptype == 0);
% assert(isempty(vectors));
% assert(AppData.S.n_base_params == Pdegree);
% assert(length(PV) == Pdegree);
% Baseline.minx = min(AppData.S.fitdata(:,5));
% maxx = max(AppData.S.fitdata(:,6));
% Baseline.XM = ((Baseline.minx:maxx)'*ones(1,Pdegree)/Pscale) .^ ...
%     (ones(maxx-Baseline.minx+1,1)*(0:Pdegree-1));
% Baseline.PVi = AppData.S.n_input_params+(1:AppData.S.n_base_params);

% function recalc = recalculate(AppData, Index)
% % recalc = recalculate(AppData, Index)
% % Produces an output matrix similar to the verbose fit output
% % matrix.
% assert(Index <= length(AppData.Scans));
% ScanNum = AppData.Scans(Index);
% IScan = AppData.IScans(Index);
% PScan = AppData.PTEi(Index);
% S = AppData.S;
% assert(ScanNum == S.scannum(IScan));
% ifile = mlf_path(AppData.ScanBase,ScanNum,'.dat');
% raw = loadbin(ifile); % Includes zero
% raw = raw - mean(raw(AppData.BGi,1));
% SR = (S.fitdata(IScan,5):S.fitdata(IScan,6))';
% BP = S.fitdata(IScan,AppData.Baseline.PVi)';
% baseline = AppData.Baseline.XM(SR - AppData.Baseline.minx+1,:) * BP;
% nu_rel = -S.EtalonFSR*etln_evalJ(AppData.PTE(PScan,5:11), ...
%     (SR-AppData.PTE(PScan,4)+1)/AppData.Baseline.Pscale);
% Abs = zeros(size(SR));
% XK = zeros(length(SR),2*S.n_lines);
% for i = 1:S.n_lines
%     X = (nu_rel + S.nu_F0(IScan) - S.nu_P(IScan,i)) / S.Ged(IScan,i);
%     Y = S.Gl(IScan,i)/S.Ged(IScan,i);
%     [K,~,~] = humdev(X,Y);
%     XK(:,2*i-1) = X;
%     XK(:,2*i) = K;
%     Abs = Abs + S.Nfit(IScan,i)*S.Scorr(IScan,i)*S.CavLen*K / ...
%         (S.Ged(IScan,i) * sqrt(pi));
% end
% fit = baseline .* exp(-S.N_Passes * Abs);
% recalc = [ SR nu_rel raw(SR,1) fit baseline Abs XK ];

function [fit, recalc] = load_fit_and_recalc(AppData, Index)
% Load verbose fit and recalc. Aborts if points don't match
recalc = icosfit_recalculate(AppData.RC, Index);
ScanNum = AppData.Scans(Index);
if AppData.S.ICOS_debug
    ffile = sprintf( '%s/%04d.dat', AppData.S.base, ScanNum);
else
    ffile = mlf_path( AppData.S.base, ScanNum, '.dat');
end
fit = load(ffile);
assert(size(recalc,1) == size(fit,1));
assert(all(recalc(:,1) == fit(:,1)));
nfc = size(fit,2);
if nfc < size(recalc,2)
    recalc = recalc(:,1:nfc);
end
assert(size(recalc,2) == size(fit,2));

function [maxabsdev, stddev, maxval] = blind_test(AppData)
maxabsdev = zeros(length(AppData.Scans),6);
stddev = zeros(size(maxabsdev));
maxval = stddev;
for i = 1:length(AppData.Scans)
    [fit, recalc] = load_fit_and_recalc(AppData, i);
    dev = fit - recalc;
    maxabsdev(i,:) = max(abs(dev));
    stddev(i,:) = std(dev);
    maxval(i,:) = max(fit);
end
% maxmaxabsdev = max(maxabsdev)
% maxstddev = max(stddev)


function recalcview_callback(handles, sv_axes)
if nargin < 2
    sv_axes = handles.Axes;
end
AppData = handles.data.AppData;
scan = handles.data.Scans(handles.data.Index);
[fe, recalc] = load_fit_and_recalc(AppData, handles.data.Index);
data_ok = (~isempty(fe));
if ~isfield(AppData,'menus')
    top_menu = uimenu(handles.scan_viewer,'Tag','recalc','Label','recalc');
    cb = @recalc_menu_callback;
    AppData.menus.Y_nu = ...
        uimenu(top_menu,'Tag','Y_nu','Label','nu','Callback',cb);
    AppData.menus.Y_raw = ...
        uimenu(top_menu,'Tag','Y_raw','Label','Raw','Callback',cb);
    AppData.menus.Y_fit = ...
        uimenu(top_menu,'Tag','Y_fit','Label','Fit','Checked','on','Callback',cb);
    AppData.menus.Y_baseline = ...
        uimenu(top_menu,'Tag','Y_baseline','Label','Baseline','Callback',cb);
    AppData.menus.Y_abs = ...
        uimenu(top_menu,'Tag','Y_abs','Label','Abs','Callback',cb);
    if size(fe,2) > 6
        AppData.menus.XK = zeros(size(fe,2)-6,1);
        for i = 7:2:size(fe,2)-1
            tag = sprintf('YX%d', i);
            AppData.menus.XK(i-6) = ...
                uimenu(top_menu,'Tag', tag, 'Label', tag(2:end), 'Callback',cb);
            tag = sprintf('YK%d', i+1);
            AppData.menus.XK(i-6+1) = ...
                uimenu(top_menu,'Tag', tag, 'Label', tag(2:end), 'Callback',cb);
        end
    else
        AppData.menus.XK = [];
    end
    AppData.menus.Detrend = ...
        uimenu(top_menu,'Tag','Detrend','Label','Detrend','Callback',cb);
    handles.data.AppData = AppData;
    guidata(handles.scan_viewer,handles);
end
if data_ok
    % AppData.plotbase = strcmp(get(AppData.menus.Baselines,'Checked'),'on');
    S = AppData.S;
    i = find(S.scannum==scan);
    % lpos = S.nu + S.delta*(S.P_vec(i)/760.) ...
    %           - S.fitdata(i,S.v);
    % if S.nu0 ~= 0
    %   lpos = lpos + S.nu_F0(i);
    % end
    % lgd = S.fitdata(i,S.v+1);
    % lgl = S.fitdata(i,S.v+3);
    % lwid = (lgd+lgl); % fitdata(i,v+1)+fitdata(i,v+3);
    nux = fe(:,2) + S.nu_F0(i) + S.nu0;
    xdir = 'reverse';
    ttlx = 'Wavenumber (cm-1)';
    dataX = nux;
    resX = nux;
    dataYa = fe(:,AppData.Yopt);
    dataYb = recalc(:,AppData.Yopt);
    resY = dataYa - dataYb;
    sdev = sqrt(mean((resY).^2));
    reslbl = sprintf('%s Res', AppData.ylbl);
    ttltext = sprintf('Scan: %d \\sigma = %f', ...
        S.scannum(i), sdev);
   
    % X = [ lpos lpos-lwid; lpos lpos+lwid ];
    % % Y = [ basep meanp; fitp meanp ];
    % Y = [ maxp meanp; fitp meanp ];
      
    plot(sv_axes(1),resX, resY);
    set(sv_axes(1),'XDir',xdir,'xticklabel',[],'Xgrid','on','Ygrid','on');
    ylabel(sv_axes(1),reslbl);
    title(sv_axes(1),ttltext);

    if AppData.Detrend
        rYa = detrend(dataYa);
        rYb = dataYb + rYa - dataYa;
        plot(sv_axes(2), dataX, rYa, dataX, rYb);
    else
        plot(sv_axes(2), dataX, dataYa, dataX, dataYb);
    end
    set(sv_axes(2),'XDir', xdir,'YAxisLocation','right','Xgrid','on','Ygrid','on');
    xlabel(sv_axes(2),ttlx);
end

function recalc_menu_callback(hObject,~)
handles = guidata(hObject);
AppData = handles.data.AppData;
Tag = get(hObject,'Tag');
if Tag(1) == 'Y'
    nu = 'off';
    raw = 'off';
    fit = 'off';
    baseline = 'off';
    abs = 'off';
    switch hObject
        case AppData.menus.Y_nu
            nu = 'on';
            handles.data.AppData.Yopt = 2;
            handles.data.AppData.ylbl = 'nu';
        case AppData.menus.Y_raw
            raw = 'on';
            handles.data.AppData.Yopt = 3;
            handles.data.AppData.ylbl = 'raw';
        case AppData.menus.Y_fit
            fit = 'on';
            handles.data.AppData.Yopt = 4;
            handles.data.AppData.ylbl = 'fit';
        case AppData.menus.Y_baseline
            baseline = 'on';
            handles.data.AppData.Yopt = 5;
            handles.data.AppData.ylbl = 'baseline';
        case AppData.menus.Y_abs
            abs = 'on';
            handles.data.AppData.Yopt = 6;
            handles.data.AppData.ylbl = 'abs';
        otherwise
            handles.data.AppData.Yopt = str2double(Tag(3:end));
            handles.data.AppData.ylbl = Tag(2:end);
    end
    set(AppData.menus.Y_nu,'checked',nu);
    set(AppData.menus.Y_raw,'checked',raw);
    set(AppData.menus.Y_fit,'checked',fit);
    set(AppData.menus.Y_baseline,'checked',baseline);
    set(AppData.menus.Y_abs,'checked',abs);
    for i=1:length(AppData.menus.XK)
        checked = 'off';
        if handles.data.AppData.Yopt == i+6
            checked = 'on';
        end
        set(AppData.menus.XK(i),'checked',checked);
    end
elseif hObject == AppData.menus.Detrend
    handles.data.AppData.Detrend = ~AppData.Detrend;
    if handles.data.AppData.Detrend
        set(AppData.menus.Detrend,'checked','on');
    else
        set(AppData.menus.Detrend,'checked','off');
    end
end
handles.data.ylim{2} = [];
guidata(hObject,handles);
scan_viewer('scan_display',handles);
