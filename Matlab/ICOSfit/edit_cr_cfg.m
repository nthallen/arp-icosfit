function varargout = edit_cr_cfg(varargin)
% EDIT_CR_CFG M-file for edit_cr_cfg.fig
%      EDIT_CR_CFG, by itself, creates a new EDIT_CR_CFG or raises the existing
%      singleton*.
%
%      H = EDIT_CR_CFG returns the handle to a new EDIT_CR_CFG or the handle to
%      the existing singleton*.
%
%      EDIT_CR_CFG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDIT_CR_CFG.M with the given input arguments.
%
%      EDIT_CR_CFG('Property','Value',...) creates a new EDIT_CR_CFG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before edit_cr_cfg_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to edit_cr_cfg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help edit_cr_cfg

% Last Modified by GUIDE v2.5 16-Jul-2007 13:56:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @edit_cr_cfg_OpeningFcn, ...
                   'gui_OutputFcn',  @edit_cr_cfg_OutputFcn, ...
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


% --- Executes just before edit_cr_cfg is made visible.
function edit_cr_cfg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to edit_cr_cfg (see VARARGIN)

% Choose default command line output for edit_cr_cfg
if nargin >= 2 && strcmp(varargin{1},'Config')
  cr_cfg = varargin{2};
else
  cr_cfg = load_cr_cfg(1);
end
set(handles.HomeDir,'String', cr_cfg.HomeDir);
set(handles.Matlab_CD_Path,'String',cr_cfg.Matlab_CD_Path);
set(handles.ICOSfit_CD_Path,'String',cr_cfg.ICOSfit_CD_Path);
set(handles.CPCI14link,'String',cr_cfg.CPCI14link);
S{1} = pwd;
cd ..
S{2} = pwd;
cd ..
S{3} = pwd;
cd(S{1});
set(handles.SaveDir, 'String', S, 'value', 1);
handles.output = cr_cfg;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes edit_cr_cfg wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = edit_cr_cfg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);

function Matlab_CD_Path_Callback(hObject, eventdata, handles)
% hObject    handle to Matlab_CD_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Matlab_CD_Path as text
%        str2double(get(hObject,'String')) returns contents of Matlab_CD_Path as a double


% --- Executes during object creation, after setting all properties.
function Matlab_CD_Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Matlab_CD_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ICOSfit_CD_Path_Callback(hObject, eventdata, handles)
% hObject    handle to ICOSfit_CD_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ICOSfit_CD_Path as text
%        str2double(get(hObject,'String')) returns contents of ICOSfit_CD_Path as a double


% --- Executes during object creation, after setting all properties.
function ICOSfit_CD_Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ICOSfit_CD_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HomeDir_Callback(hObject, eventdata, handles)
% hObject    handle to HomeDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HomeDir as text
%        str2double(get(hObject,'String')) returns contents of HomeDir as a double


% --- Executes during object creation, after setting all properties.
function HomeDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HomeDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CPCI14link_Callback(hObject, eventdata, handles)
% hObject    handle to CPCI14link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CPCI14link as text
%        str2double(get(hObject,'String')) returns contents of CPCI14link as a double


% --- Executes during object creation, after setting all properties.
function CPCI14link_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPCI14link (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_btn.
function Save_btn_Callback(hObject, eventdata, handles)
% hObject    handle to Save_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cr_cfg.HomeDir = get(handles.HomeDir,'String');
cr_cfg.Matlab_CD_Path = get(handles.Matlab_CD_Path,'String');
cr_cfg.ICOSfit_CD_Path = get(handles.ICOSfit_CD_Path,'String');
cr_cfg.CPCI14link = get(handles.CPCI14link,'String');
handles.output = cr_cfg;
guidata(hObject, handles);
SaveDirs = get(handles.SaveDir,'String');
SaveDir = SaveDirs{get(handles.SaveDir,'value')};
fd = fopen([ SaveDir '/CR_Config.m'], 'w');
fprintf(fd, 'function cr_cfg = CR_Config;\n');
fprintf(fd, '% CR_Config defines local configuration\n');
fprintf(fd, 'cr_cfg.Matlab_CD_Path = ''%s'';\n', cr_cfg.Matlab_CD_Path );
fprintf(fd, 'cr_cfg.ICOSfit_CD_Path = ''%s'';\n', cr_cfg.ICOSfit_CD_Path );
fprintf(fd, 'cr_cfg.Homedir = ''%s'';\n', cr_cfg.HomeDir );
fprintf(fd, 'cr_cfg.CPCI14link = ''%s'';\n', cr_cfg.CPCI14link );
fclose(fd);
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);

