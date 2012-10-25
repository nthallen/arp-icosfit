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
    60    45    60     1    20    30     0     1
    60    45    60     1     0    30    60     1
    ];
scan_viewer('Scans', scans, 'Axes', Axes, 'Name', 'Scan View', ...
    'Callback', @scanview_callback, 'AppData', AppData);

function scanview_callback(handles)
scan = handles.data.Scans(handles.data.Index);
path = mlf_path( handles.data.AppData.base, scan, '.dat');
fe = loadbin( path );
data_ok = (length(fe)>0);
if data_ok
    nsamples = size(fe,1);
    plot(handles.Axes(1),[1:nsamples],fe(:,1));
    set(handles.Axes(1),'xticklabel',[]);
    title(handles.Axes(1),sprintf('Scan %d: %s', scan, getrun ));

    plot(handles.Axes(2), [1:nsamples], fe(:,2));
    set(handles.Axes(2),'YAxisLocation','right');
    xlabel(handles.Axes(2),'Samples');
end
