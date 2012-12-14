function scanview( scans, base );
% scanview( scans [, base] );
% Display raw scan with etalon, optionally limited to a range
% of Scan values.
if nargin < 1
  base = '';
  scans = listscans();
elseif nargin == 1
  base = '';
end
AppData.base = find_scans_dir(base);

Axes = [
    60    45    60     1    20    30     0     1    0
    60    45    60     1     0    30    60     1    0
    ];
scan_viewer('Scans', scans, 'Axes', Axes, 'Name', 'Scan View', ...
    'Callback', @scanview_callback, 'AppData', AppData);

function scanview_callback(handles, sv_axes)
if nargin < 2
    sv_axes = handles.Axes;
end
scan = handles.data.Scans(handles.data.Index);
path = mlf_path( handles.data.AppData.base, scan, '.dat');
fe = loadbin( path );
data_ok = (~isempty(fe));
if data_ok
    nsamples = size(fe,1);
    plot(sv_axes(1),[1:nsamples],fe(:,1));
    set(sv_axes(1),'xticklabel',[]);
    title(sv_axes(1),sprintf('Scan %d: %s', scan, ...
        getrun(0,handles.scan_viewer)));

    plot(sv_axes(2), [1:nsamples], fe(:,2));
    set(sv_axes(2),'YAxisLocation','right');
    xlabel(sv_axes(2),'Samples');
end
