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

% Last Modified by GUIDE v2.5 22-Oct-2012 16:22:34

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
handles.data.scan_increment = 0;
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


% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function ViewerGroup_SelectionChangeFcn(hObject, eventdata)
handles = guidata(hObject);
if eventdata.NewValue == handles.Play
    handles.data.scan_increment = 1;
elseif eventdata.NewValue == handles.FastFwd
    if handles.data.scan_increment > 0
        handles.data.scan_increment = handles.data.scan_increment * 3;
    else
        handles.data.scan_increment = 10;
    end
elseif eventdata.NewValue == handles.Pause
    handles.data.scan_increment = 0;
elseif eventdata.NewValue == handles.Reverse
    handles.data.scan_increment = -1;
elseif eventdata.NewValue == handles.FastRev
    if handles.data.scan_increment < 0
        handles.data.scan_increment = handles.data.scan_increment * 3;
    else
        handles.data.scan_increment = -10;
    end
else
    errordlg('Unknown object');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CrntScan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrntScan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
