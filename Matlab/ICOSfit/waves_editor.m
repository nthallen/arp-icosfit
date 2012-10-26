function varargout = waves_editor(varargin)
% WAVES_EDITOR M-file for waves_editor.fig
% WAVES_EDITOR( 'data', data, 'index', 3, 'qclicomp', 1 );
% WAVES_EDITOR( 'cpci', cpci_number[, 'qclicomp', 1] );
%      data is a vector of raw data
%      index specifies which waveform to start editting
%      qclicomp if specified limits which elements of the waveform can be
%        edited, and does not allow switching between waveforms, since the
%        input data is presumably from a specific pre-defined waveform.
%      If cpci is specified, the data will be retrieved from the usual
%        place and the waveform will be identified from PT.mat. In this
%        case, qclicomp is assumed to be true.

%      WAVES_EDITOR, by itself, creates a new WAVES_EDITOR or raises the existing
%      singleton*.
%
%      H = WAVES_EDITOR returns the handle to a new WAVES_EDITOR or the handle to
%      the existing singleton*.
%
%      WAVES_EDITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVES_EDITOR.M with the given input arguments.
%
%      WAVES_EDITOR('Property','Value',...) creates a new WAVES_EDITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before waves_editor_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to waves_editor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% TODO:
%  Delete (X) does not close
%  Provide a read-only waves.m version that only updates SignalRegion
%    would disable most selections, including TriggerDelay
%  Move axes down a bit so exponents show up

% Edit the above text to modify the response to help waves_editor

% Last Modified by GUIDE v2.5 13-Feb-2007 11:35:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @waves_editor_OpeningFcn, ...
                   'gui_OutputFcn',  @waves_editor_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before waves_editor is made visible.
function waves_editor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to waves_editor (see VARARGIN)

% Choose default command line output for waves_editor
handles.output = hObject;

% Check for 'data' property
handles.data.index = 1;
handles.data.rawdata = [];
handles.data.qclicomp = 0;
for i=1:2:length(varargin)-1
  if strcmp(varargin{i},'data')
    handles.data.rawdata = varargin{i+1};
  elseif strcmp(varargin{i},'index')
    handles.data.index = varargin{i+1};
  elseif strcmp(varargin{i},'qclicomp')
    handles.data.qclicomp = varargin{i+1};
  elseif strcmp(varargin{i}, 'cpci')
    handles.data.qclicomp = 1;
    cpci = varargin{i+1};
    handles.data.rawdata = scanload(cpci);
    PT = load_mat_files('PT');
    v = find(diff(PT.CPCI14)>0)+1;
    wvno = interp1(PT.CPCI14(v), PT.QCLI_Wave(v), cpci, 'nearest')+1;
    handles.data.index = wvno;
  end
end

% Load waves.m if it exists, and initialize handles and waves_menu
handles.data.WaveSpecs = load_waves;
% if exist(ICOSfit,'file') == 2
%   waves;
%   handles.data.WaveSpecs = load_waves;
% else
%   if handles.data.qclicomp
%     error('Unable to locate waves.m');
%   end
%   handles.data.WaveSpecs = struct( ...
%     'Name', { 'Basic' }, ...
%     'RawRate', { 100000 }, ...
%     'RawSamples', { 1000 }, ...
%     'NetSamples', { 1000 }, ...
%     'NAverage', { 1 }, ...
%     'NCoadd', { 10 }, ...
%     'FTrigger', { 100 }, ...
%     'ISICOS', { 1 }, ...
%     'TzSamples', { 350 }, ...
%     'dIdt', { 30 }, ...
%     'I0', { 0.64 } ...
%     );
%   handles.data.index = 1;
% end
handles.data.xwvs = {};
if handles.data.qclicomp
  bgcolor = get(handles.uipanel1,'BackgroundColor');
  set(handles.waves_menu,'enable','off');
  set(handles.RawRate,'enable','off','BackgroundColor',bgcolor);
  set(handles.New_Wave_Btn,'Enable','off','visible','off');
  set(handles.NAverage,'enable','off','BackgroundColor',bgcolor);
  set(handles.NetSamples,'enable','off','BackgroundColor',bgcolor);
  set(handles.TriggerDelay,'Enable','off','visible','off');
  set(handles.TDlabel','visible','off');
  set(handles.NCoadd,'enable','off','BackgroundColor',bgcolor);
  set(handles.FTrigger,'enable','off','BackgroundColor',bgcolor);
  set(handles.TzSamples,'enable','off','BackgroundColor',bgcolor);
end

% Update handles structure
guidata(hObject, handles);

% Call the function that updates all of the parameters on the screen
setup_waveform(hObject, handles );


% UIWAIT makes waves_editor wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = waves_editor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.data.index;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes when user attempts to close figure1.
function waves_editor_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end


% --- Executes on selection change in waves_menu.
function waves_menu_Callback(hObject, eventdata, handles)
% hObject    handle to waves_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns waves_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from waves_menu
index = get(hObject,'Value');
handles.data.index = index;
guidata(hObject,handles);
setup_waveform(hObject, handles);

% --- Executes during object creation, after setting all properties.
function waves_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to waves_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in New_Wave_Btn.
function New_Wave_Btn_Callback(hObject, eventdata, handles)
% hObject    handle to New_Wave_Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newname = inputdlg('New Waveform Name:', 'New Waveform');
if length(newname) && length(newname{1})
  if any(strcmp(newname{1},{handles.data.WaveSpecs.Name}))
    errordlg(sprintf('A waveform named %s already exists', newname{1}));
  else
    index = handles.data.index;
    wv = handles.data.WaveSpecs(index);
    if index <= length(handles.data.xwvs)
      xwv = handles.data.xwvs{handles.data.index};
    else
      xwv = [];
    end
    wv.Name = newname{1};
    handles.data.WaveSpecs(end+1) = wv;
    handles.data.index = length(handles.data.WaveSpecs);
    handles.data.xwvs{handles.data.index} = xwv;
    guidata(hObject, handles);
    setup_waveform(hObject, handles);
  end
end


function RawRate_Callback(hObject, eventdata, handles)
% hObject    handle to RawRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RawRate as text
%        str2double(get(hObject,'String')) returns contents of RawRate as a double


% --- Executes during object creation, after setting all properties.
function RawRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RawRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RawSamples_Callback(hObject, eventdata, handles)
% hObject    handle to RawSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RawSamples as text
%        str2double(get(hObject,'String')) returns contents of RawSamples as a double


% --- Executes during object creation, after setting all properties.
function RawSamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RawSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NAverage_Callback(hObject, eventdata, handles)
% hObject    handle to NAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NAverage as text
%        str2double(get(hObject,'String')) returns contents of NAverage as a double


% --- Executes during object creation, after setting all properties.
function NAverage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NCoadd_Callback(hObject, eventdata, handles)
% hObject    handle to NCoadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NCoadd as text
%        str2double(get(hObject,'String')) returns contents of NCoadd as a double


% --- Executes during object creation, after setting all properties.
function NCoadd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NCoadd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to FTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FTrigger as text
%        str2double(get(hObject,'String')) returns contents of FTrigger as a double


% --- Executes during object creation, after setting all properties.
function FTrigger_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TzSamples_Callback(hObject, eventdata, handles)
% hObject    handle to TzSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TzSamples as text
%        str2double(get(hObject,'String')) returns contents of TzSamples as a double


% --- Executes during object creation, after setting all properties.
function TzSamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TzSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NetSamples_Callback(hObject, eventdata, handles)
% hObject    handle to NetSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NetSamples as text
%        str2double(get(hObject,'String')) returns contents of NetSamples as a double


% --- Executes during object creation, after setting all properties.
function NetSamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NetSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_Btn.
function Save_Btn_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
WaveSpecs = handles.data.WaveSpecs;
if ~handles.data.qclicomp
  save waves.mat WaveSpecs;
  fid = fopen( 'waves.m', 'w' );
  fprintf( fid, 'load(''waves.mat'');\n' );
  fclose(fid);
end
% Still need to do the <wave>_etln.mat file
for index=1:length(WaveSpecs)
  if index <= length(handles.data.xwvs) && length(handles.data.xwvs{index})
    update_etln_file(WaveSpecs(index).Name, handles.data.xwvs{index});
  end
end
uiresume(handles.figure1);

% --- Executes on button press in Cancel_Btn.
function Cancel_Btn_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


function SignalStart_Callback(hObject, eventdata, handles)
% hObject    handle to SignalStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SignalStart as text
%        str2double(get(hObject,'String')) returns contents of SignalStart as a double


% --- Executes during object creation, after setting all properties.
function SignalStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SignalStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SignalEnd_Callback(hObject, eventdata, handles)
% hObject    handle to SignalEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SignalEnd as text
%        str2double(get(hObject,'String')) returns contents of SignalEnd as a double


% --- Executes during object creation, after setting all properties.
function SignalEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SignalEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Pick_Tz.
function Pick_Tz_Callback(hObject, eventdata, handles)
% hObject    handle to Pick_Tz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Pick_SR.
function Pick_SR_Callback(hObject, eventdata, handles)
% hObject    handle to Pick_SR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function NScans_Callback(hObject, eventdata, handles)
% hObject    handle to NScans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NScans as text
%        str2double(get(hObject,'String')) returns contents of NScans as a double


% --- Executes during object creation, after setting all properties.
function NScans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NScans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InputSamples_Callback(hObject, eventdata, handles)
% hObject    handle to InputSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputSamples as text
%        str2double(get(hObject,'String')) returns contents of InputSamples as a double


% --- Executes during object creation, after setting all properties.
function InputSamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function setup_waveform(hObject, handles )
index = handles.data.index;
set(handles.waves_menu,'String',{handles.data.WaveSpecs.Name},'value',index);
wv = handles.data.WaveSpecs(index);
if index <= length(handles.data.xwvs) && length(handles.data.xwvs{index})
  xwv = handles.data.xwvs{handles.data.index};
else
  xwv = get_waveform_params( wv.Name, 'TriggerDelay', 0, ...
    'holdoff', 4e-4, 'CurRamp', 1 );
  %   fname = findinpath( [ wv.Name '_etln.mat' ], { '.', '..', '../..' } );
  %   if length(fname)
  %     xwv = load(fname);
  %   else
  %     xwv = struct('TriggerDelay',0);
  %   end
end
if isfield(xwv,'SignalRegion')
  startsample = min(xwv.SignalRegion);
  endsample = max(xwv.SignalRegion);
else
  % xwv.holdoff = 4e-4;
  startsample = round(wv.TzSamples + xwv.holdoff*round(wv.RawRate/wv.NAverage));
  endsample = wv.NetSamples - wv.TzSamples - 1;
  xwv.SignalRegion = startsample:endsample;
end
% if ~isfield(xwv,'TriggerDelay')
%   xwv.TriggerDelay = 0;
% end
% if isfield(xwv,'RampRegions' ) && isfield(xwv, 'CurRamp' )
if isfield(xwv,'RampRegions' )
  n_ramps = size(xwv.RampRegions,1);
  if xwv.CurRamp < 1 || xwv.CurRamp > n_ramps
    xwv.CurRamp = 1;
  end
else
  n_ramps = 1;
  xwv.CurRamp = 1;
end
set(handles.CurRamp,'String',cellstr(num2str((1:n_ramps)')), ...
  'value',xwv.CurRamp);
if n_ramps > 1
  set(handles.CurRamp, 'enable', 'on');
else
  set(handles.CurRamp, 'enable', 'off');
end

handles.data.xwvs{index} = xwv;
guidata(hObject, handles);
nraw = length(handles.data.rawdata);
NScans = floor(nraw/wv.NetSamples);
if nraw
  set(handles.InputSamples,'String',num2str(nraw));
  if NScans*wv.NetSamples == nraw
    set(handles.NScans,'String',num2str(NScans),'backgroundcolor', ...
      get(handles.uipanel1,'backgroundcolor'));
  else
    set(handles.NScans,'String',num2str(nraw/wv.NetSamples), ...
      'backgroundcolor',[1 1 .5]);
  end
else
  set(handles.InputSamples,'String',num2str(wv.NetSamples));
  set(handles.NScans,'String','1');
  nraw = 0;
end
set(handles.index,'String',num2str(handles.data.index-1));
set(handles.RawRate,'String',num2str(wv.RawRate));
set(handles.RawSamples,'String',num2str(wv.RawSamples));
set(handles.NAverage,'String',num2str(wv.NAverage));
set(handles.NetSamples,'String',num2str(wv.NetSamples));
set(handles.NCoadd,'String',num2str(wv.NCoadd));
set(handles.FTrigger,'String',num2str(wv.FTrigger));
set(handles.TzSamples,'String',num2str(wv.TzSamples));
set(handles.TriggerDelay,'String',num2str(xwv.TriggerDelay));
set(handles.SignalStart,'String',num2str(startsample));
set(handles.SignalEnd,'String',num2str(endsample));
if nraw
  D = handles.data.rawdata;
  if handles.data.qclicomp == 0 && xwv.TriggerDelay
    D = [ D(xwv.TriggerDelay+1:end); D(1:xwv.TriggerDelay) ];
  end
  if NScans > 1
    E = reshape(D(1:wv.NetSamples*NScans),wv.NetSamples,NScans);
  else
    E = D;
  end
  axes(handles.axes1);
  plot(E);
  yl = ylim;
  tz = wv.TzSamples;
  ys = [yl(1) yl(2) yl(2) yl(1)];
  h = fill([1 1 tz tz],ys,[.8 .8 1]);
  % set(h,'LineStyle','none');
  hold on;
  h = fill([startsample startsample endsample endsample],ys,[.8 1 .8]);
  % set(h,'LineStyle','none');
  plot(E);
  hold off;
  %set(handles.axes1,'layer','top');
end

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Apply_Btn.
function Apply_Btn_Callback(hObject, eventdata, handles)
% hObject    handle to Apply_Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index = str2num(get(handles.index,'String'))+1;
wv = handles.data.WaveSpecs(index);
wv.RawRate = str2num(get(handles.RawRate,'String'));
wv.NAverage = str2num(get(handles.NAverage,'String'));
wv.NetSamples = str2num(get(handles.NetSamples,'String'));
wv.RawSamples = wv.NetSamples*wv.NAverage;
set(handles.RawSamples,'String',num2str(wv.RawSamples));
wv.NCoadd = str2num(get(handles.NCoadd,'String'));
wv.FTrigger = str2num(get(handles.FTrigger,'String'));
wv.TzSamples = str2num(get(handles.TzSamples,'String'));
handles.data.WaveSpecs(index) = wv;
xwv = handles.data.xwvs{index};
xwv.TriggerDelay = str2num(get(handles.TriggerDelay,'String'));
startsample = str2num(get(handles.SignalStart,'String'));
endsample = str2num(get(handles.SignalEnd,'String'));
xwv.SignalRegion = startsample:endsample;
if isfield(xwv, 'RampRegions') && isfield(xwv, 'CurRamp')
  xwv.RampRegions(xwv.CurRamp,:) = [ startsample endsample ];
end
handles.data.xwvs{index} = xwv;
guidata(hObject, handles);
setup_waveform(hObject, handles);

function update_etln_file( waveform, xwv )
save([waveform '_etln.mat'], '-struct', 'xwv' );

function Holdoff_CreateFcn(hObject, eventdata, handles)
return;

function Holdoff_Callback(hObject, eventdata, handles)
return;



% --- Executes on button press in AddRamp.
function AddRamp_Callback(hObject, eventdata, handles)
% hObject    handle to AddRamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
index = str2num(get(handles.index,'String'))+1;
xwv = handles.data.xwvs{index};
if ~isfield(xwv, 'RampRegions')
  xwv.RampRegions = [min(xwv.SignalRegion) max(xwv.SignalRegion)];
end
xwv.RampRegions = [
  xwv.RampRegions;
  [min(xwv.SignalRegion) max(xwv.SignalRegion)] ];
xwv.CurRamp = size(xwv.RampRegions,1);
handles.data.xwvs{index} = xwv;
guidata(hObject, handles);
setup_waveform( hObject, handles );

% --- Executes on selection change in CurRamp.
function CurRamp_Callback(hObject, eventdata, handles)
% hObject    handle to CurRamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns CurRamp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CurRamp
index = str2num(get(handles.index,'String'))+1;
xwv = handles.data.xwvs{index};
n_ramps = size(xwv.RampRegions,1);
CurRamp = get(hObject,'Value');
if CurRamp < 1 || CurRamp > n_ramps
  errdialog('Invalid CurRamp value');
else
  xwv.CurRamp = CurRamp;
  xwv.SignalRegion = ...
    [ xwv.RampRegions(CurRamp,1):xwv.RampRegions(CurRamp,2) ];
  handles.data.xwvs{index} = xwv;
  guidata( hObject, handles );
  setup_waveform( hObject, handles );
end


% --- Executes during object creation, after setting all properties.
function CurRamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurRamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


