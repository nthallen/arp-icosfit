function varargout = scan_viewer(varargin)
% SCAN_VIEWER MATLAB code for scan_viewer.fig
%      SCAN_VIEWER, by itself, creates a new SCAN_VIEWER or raises the existing
%      singleton*.
%
%      H = SCAN_VIEWER returns the handle to a new SCAN_VIEWER or the handle to
%      the existing singleton*.
%
%      SCAN_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCAN_VIEWER.M with the given input arguments.
%
%      SCAN_VIEWER('Property','Value',...) creates a new SCAN_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scan_viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scan_viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scan_viewer

% Last Modified by GUIDE v2.5 22-Oct-2012 21:53:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @scan_viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @scan_viewer_OutputFcn, ...
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


% --- Executes just before scan_viewer is made visible.
function scan_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scan_viewer (see VARARGIN)

% Choose default command line output for scan_viewer
handles.output = hObject;

% Update handles structure
set(handles.ViewerGroup, 'SelectionChangeFcn', ...
    @ViewerGroup_SelectionChangeFcn);
handles.data.Scans = [];
handles.data.Index = 1;
for i=1:2:length(varargin)-1
  if strcmpi(varargin{i},'Scans')
      handles.data.Scans = varargin{i+1};
  else
      errordlg(sprintf('Unrecognized property: %s', varargin{i}));
  end
end
handles.data.Index_max = length(handles.data.Scans);
if isempty(handles.data.Scans)
    errordlg('No scans specified','Trouble','modal');
else
    set(handles.Slider,'Min',1,'Max',handles.data.Index_max,'Value',1);
    scan_display(handles);
end
guidata(hObject, handles);

% UIWAIT makes scan_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = scan_viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on Slider movement.
function Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of Slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of Slider
opmode = get(handles.ViewerGroup,'SelectedObject');
if opmode == handles.Pause
    handles.data.Index = round(get(hObject,'Value'));
    scan_display(handles);
    guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function Slider_CreateFcn(hObject, ~, ~)
% hObject    handle to Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: Slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function ViewerGroup_SelectionChangeFcn(hObject, ~)
while 1
    handles = guidata(hObject);
    Speed = get(handles.Speed,'UserData');
    action = get(handles.ViewerGroup,'SelectedObject');
    SIndex = round(get(handles.Slider,'Value'));
    if SIndex ~= handles.data.Index && action ~= handles.Pause
        set(handles.ViewerGroup,'SelectedObject',handles.Pause);
        handles.data.Index = SIndex;
        scan_display(handles);
    end
    if action == handles.Pause
        break;
    elseif action == handles.Play
        if handles.data.Index < handles.data.Index_max
            handles.data.Index = handles.data.Index+1;
            scan_display(handles);
        end
        set(handles.ViewerGroup,'SelectedObject',handles.Pause);
    elseif action == handles.FastFwd
        if handles.data.Index + Speed(2) < handles.data.Index_max
            handles.data.Index = handles.data.Index + Speed(2);
            scan_display(handles);
            pause(Speed(1));
        else
            set(handles.ViewerGroup,'SelectedObject',handles.Pause);
        end
    elseif action == handles.Reverse
        if handles.data.Index > 1
            handles.data.Index = handles.data.Index-1;
            scan_display(handles);
        end
        set(handles.ViewerGroup,'SelectedObject',handles.Pause);
    elseif action == handles.FastRev
        if handles.data.Index > 1
            if handles.data.Index > Speed(2)
                handles.data.Index = handles.data.Index - Speed(2);
            else
                handles.data.Index = 1;
            end
            scan_display(handles);
            pause(Speed(1));
        else
            set(handles.ViewerGroup,'SelectedObject',handles.Pause);
        end
    else
        errordlg('Unknown object');
    end
    guidata(hObject, handles);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CrntScan_CreateFcn(~, ~, ~)
% hObject    handle to CrntScan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function P = Set_Speed(hObject, handle, Pin, Pnew)
if hObject == handle
    set(handle,'checked','on');
    P = Pnew;
else
    set(handle,'checked','off');
    P = Pin;
end

% --------------------------------------------------------------------
function Speed_Callback(hObject, ~, handles)
% hObject    handle to Speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
P = get(handles.Speed,'UserData');
P = Set_Speed(hObject, handles.Spd_1Hz, P, [ 1, 1 ]);
P = Set_Speed(hObject, handles.Spd_3Hz, P, [ 1/3 1]);
P = Set_Speed(hObject, handles.Spd_10Hz, P, [1/10 1]);
P = Set_Speed(hObject, handles.Spd_Step_1, P, [0 1]);
P = Set_Speed(hObject, handles.Spd_Step_10, P, [0 10]);
P = Set_Speed(hObject, handles.Spd_Step_100, P, [0 100]);
set(handles.Speed, 'UserData', P);

function scan_display(handles)
% Update CrntScan, slider
if handles.data.Index_max >= 1
    set(handles.CrntScan,'String',num2str(handles.data.Scans(handles.data.Index)));
    set(handles.Slider,'Value',handles.data.Index);
    drawnow;
end
