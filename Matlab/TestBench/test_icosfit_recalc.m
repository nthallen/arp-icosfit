function test_icosfit_recalc(base, varargin)
% test_icosfit_recalc(base)
% test_icosfit_recalc(base, 'view')

% Tests:
%   Background: Verify that BackgroundRegion is 5:TzSamples
%   Raw: Verify that raw scan and raw column of verbose output
%     are identical
%   Base: Verify that recalculated baseline matches verbose baseline
%   nu: Verify that nu gets calculated correctly
%   Fit: Verify that recalculated fit matches verbose fit
AppData.S = ICOS_setup(base);
if AppData.S.ICOSfit_format_ver < 2
    error('ICOSfit recalculation requires ICOSfit_format_ver >= 2');
end
AppData.ICOSfit_cfg = load_ICOSfit_cfg;
AppData.ScanBase = find_scans_dir([]);
AppData.Scans = AppData.S.scannum;
AppData.IScans = 1:length(AppData.Scans);
% now AppData.Scans == AppData.S.scannum(AppData.IScans)
% When scans are specified explicitly, this allows us to deal with
% fits that might run in reverse or even have multiple overlapping
% regions.

% Load baseline stuff
AppData.Baseline = loadetlnbase(AppData);

% Load PTE stuff
if ~exist(AppData.S.PTEfile,'file')
    error('Unable to locate PTEfile %s', AppData.S.PTEfile);
end
AppData.PTE = load(AppData.S.PTEfile);

AppData.opt_view = 0; % Show results interactively
AppData.opt_raw = 0;  % Check/show that raw scans match
AppData.opt_base = 0; % Check/show that baselines match
AppData.opt_fit = 0;  % Check/show that fits match
i = 1;
while i <= length(varargin)
    if ~ischar(varargin{i})
        error('Expected option string at argument %d', i+1);
    elseif strcmpi(varargin{i}, 'view')
        AppData.opt_view = 1;
    elseif strcmpi(varargin{i}, 'raw')
        AppData.opt_raw = 1;
    elseif strcmpi(varargin{i}, 'base')
        AppData.opt_base = 1;
    elseif strcmpi(varargin{i}, 'fit')
        AppData.opt_fit = 1;
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

if ~(AppData.opt_raw || AppData.opt_base || AppData.opt_fit)
    AppData.opt_raw = 1;
    AppData.opt_base = 1;
    AppData.opt_fit = 1;
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
        'Axes', [ 60    45    60     1    20    30    60     1    0 ], ...
        'Name', 'Test ICOSfit Recalc Viewer', ...
        'Callback', @recalcview_callback, 'AppData', AppData);
else
    blind_test(AppData);
end

function Baseline = loadetlnbase(AppData)
% Assign Baseline params to structure
[ nu, vectors, Pdegree, Ptype, PV, Pscale ] = ...
    readetlnbase( AppData.S.BaselineFile );
Baseline = struct('nu', nu, 'vectors', vectors, ...
    'Pdegree', Pdegree, 'Ptype', Ptype, 'PV', PV, ...
    'Pscale', Pscale);
assert(Ptype == 0);
assert(isempty(vectors));
assert(AppData.S.n_base_params == Pdegree);
assert(length(PV) == Pdegree);
Baseline.minx = min(AppData.S.fitdata(:,5));
maxx = max(AppData.S.fitdata(:,6));
Baseline.XM = ((Baseline.minx:maxx)'*ones(1,Pdegree)/Pscale) .^ ...
    (ones(maxx-Baseline.minx+1,1)*(0:Pdegree-1));
Baseline.PVi = AppData.S.n_input_params+(1:AppData.S.n_base_params);

function recalc = recalculate(AppData, Index)
% recalc = recalculate(AppData, Index)
% Produces an output matrix similar to the verbose fit output
% matrix.
assert(Index <= length(AppData.Scans));
ScanNum = AppData.Scans(Index);
IScan = AppData.IScans(Index);
PScan = AppData.PTEi(Index);
assert(ScanNum == AppData.S.scannum(IScan));
ifile = mlf_path(AppData.ScanBase,ScanNum,'.dat');
raw = loadbin(ifile); % Includes zero
raw = raw - mean(raw(AppData.BGi,1));
SR = (AppData.S.fitdata(IScan,5):AppData.S.fitdata(IScan,6))';
BP = AppData.S.fitdata(IScan,AppData.Baseline.PVi)';
baseline = AppData.Baseline.XM(SR - AppData.Baseline.minx+1,:) * BP;
nu_rel = -AppData.S.EtalonFSR*etln_evalJ(AppData.PTE(PScan,5:11), ...
    (SR-AppData.PTE(PScan,4))/AppData.Baseline.Pscale);
recalc = [ SR nu_rel raw(SR,1) baseline ];

function [recalc, fit] = load_fit_and_recalc(AppData, Index)
% Load verbose fit and recalc. Aborts if points don't match
recalc = recalculate(AppData, Index);
ScanNum = AppData.Scans(Index);
ffile = mlf_path(AppData.S.base,ScanNum,'.dat');
fit = load(ffile);
assert(size(recalc,1) == size(fit,1));
assert(all(recalc(:,1) == fit(:,1)));

function blind_test(AppData)
maxabsdev = zeros(length(AppData.Scans),4);
stddev = zeros(size(maxabsdev));
for i = 1:length(AppData.Scans)
    [recalc, fit] = load_fit_and_recalc(AppData, i);
    dev = fit(:,[1 2 3 5]) - recalc(:,[1 2 3 4]);
    maxabsdev(i,:) = max(abs(dev));
    stddev(i,:) = std(dev);
end
maxmaxabsdev = max(maxabsdev)
maxstddev = max(stddev)


function recalcview_callback(handles, sv_axes)
if nargin < 2
    sv_axes = handles.Axes;
end
AppData = handles.data.AppData;
if ~isfield(AppData,'menus')
end
scan = handles.data.Scans(handles.data.Index);
if AppData.S.ICOS_debug
    path = sprintf( '%s/%04d.dat', AppData.S.base, scan);
else
    path = mlf_path( AppData.S.base, scan, '.dat');
end
fe = load( path );
data_ok = (~isempty(fe));
if data_ok
end

