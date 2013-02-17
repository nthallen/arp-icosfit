function varargout = etln_fit(varargin)
% ETLN_FIT M-file for etln_fit.fig
% ETLN_FIT( 'SCANNUM', scannum_vector, 'OFILE', output_filename, 'SAVEALL', 1 );
%      ETLN_FIT, by itself, creates a new ETLN_FIT or raises the existing
%      singleton*.
%
%      H = ETLN_FIT returns the handle to a new ETLN_FIT or the handle to
%      the existing singleton*.
%
%      ETLN_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ETLN_FIT.M with the given input arguments.
%
%      ETLN_FIT('Property','Value',...) creates a new ETLN_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before etln_fit_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to etln_fit_OpeningFcn via varargin.
%
%      SCANNUM - vector of scan numbers to fit
%      OFILE - output file name
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% TODO:
%   Write Output
%   Provide 'Save Defaults' to save X and possibly threshold
%   Provide a progress meter, and/or quality of fit trend

%### Trouble in 070625.2 at 337 and 852

% Last Modified by GUIDE v2.5 02-Nov-2012 15:36:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @etln_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @etln_fit_OutputFcn, ...
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


% --- Executes just before etln_fit is made visible.
function etln_fit_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to etln_fit (see VARARGIN)

% Choose default command line output for etln_fit
handles.output = hObject;
handles.data.ofile = 'PTE.txt';
handles.data.saveall = 0;
for i=1:2:length(varargin)-1
  if strcmpi(varargin{i},'SCANNUM')
    handles.data.scans = varargin{i+1};
  elseif strcmpi(varargin{i},'OFILE')
    handles.data.ofile = varargin{i+1};
  elseif strcmpi(varargin{i},'SAVEALL')
    handles.data.saveall = varargin{i+1};
  end
end
wv = waves_used(handles.data.scans);
if length(wv) > 1
  errordlg( 'Multiple waveforms specified' );
  delete(handles.figure1);
  return;
end
handles.data.wv = wv;
handles.data.vZ = 1:wv.TzSamples;
% Pick out P and T for output
PT = load_mat_files('PT');
dcpi = find(diff(PT.ScanNum)>0)+1; % index of new scannum numbers
dcpi = [ dcpi(1)-1; dcpi ];
idx = ceil(interp1(PT.ScanNum(dcpi),dcpi,handles.data.scans));
if any(isnan(idx))
  errordlg('Input scan number range exceeds PT record');
  delete(handles.figure1);
  return;
end
handles.data.P = PT.CellP(idx);
handles.data.T = PT.Tavg(idx);

handles.data.ofd = fopen( handles.data.ofile, 'a' );
if handles.data.ofd == -1
  errordlg(sprintf('Cannot open output file: %s', handles.data.ofile));
  delete(handles.figure1);
  return;
end

% Get existing configuration
[ prefilterwidth, X, range_dflt, threshold ] = get_waveform_params( wv.Name, ...
  'prefilterwidth', 5, ...
  'X', [10.2817    48.3    0   -2.6921  .1645924   -3.7796   .0689779 ], ...
  'SignalRegion', wv.TzSamples+100:wv.NetSamples-wv.TzSamples-20, ...
  'threshold', .07  );
handles.data.prefilterwidth = prefilterwidth;
set(handles.prefilterwidth,'String',num2str(prefilterwidth));
handles.data.ScanNumdir = find_scans_dir([]);
handles.data.indexes = 1:length(handles.data.scans);
handles.data.peakx = [];
handles.data.Xdflt = X;
handles.data.Xlast = [];
handles.data.Ylast = [];
handles.data.X = X;
handles.data.Y = [];
handles.data.samples = range_dflt;
handles.data.threshold = threshold;
set(handles.threshold,'String',num2str(handles.data.threshold));
handles.data.figerrs = nan * handles.data.scans;
handles.data.passes = 0*handles.data.scans;
handles.data.index = 0;
handles.data.level = 0;
handles.data.startlevel = 1;
handles.data.figerr = -1;
handles.data.stop_on_fail = 1;
handles.data.peakpts = [];
handles.data.rxs = (1:length(range_dflt))'*1e-3;
handles.data.Op = optimset('lsqcurvefit');
handles.data.Op = ...
  optimset(handles.data.Op,'Jacobian', 'on','TolFun',.1, ...
  'MaxFunEvals',100,'Display','off','Algorithm',{'levenberg-marquardt',1});

if handles.data.X(6) == 0
    set(handles.dblexp,'Value',0);
else
    set(handles.dblexp,'Value',1);
end
handles.data.polling = 0;
handles.data.Xcolor.full = [1 1 1];
handles.data.Xcolor.fit = [.8 .8 .8];
handles.data.Xcolor.fixed = [.7 .7 .7];
set(handles.tc_select,'SelectedObject',handles.X5);
% handles.data.SliderSelected = handles.X5;
handles.SliderListener = ...
    addlistener(handles.Slider,'Value','PostSet', ...
    @(hObj,eventdat)etln_fit('Slider_Motion', ...
       hObject,eventdat,guidata(hObject)));
warning('off','MATLAB:rankDeficientMatrix');
set(handles.Fitting,'visible','off');
% Update handles structure
guidata(hObject, handles);
update_X_to_fig(handles);
dblexp_Callback(handles.dblexp, [], handles);

% handles = setup_level(hObject, handles);
next_scannum_file(hObject, handles );


% UIWAIT makes etln_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = etln_fit_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% --- Executes on button press in handpickpeaks_btn.
function handpickpeaks_btn_Callback(hObject, ~, handles)
% hObject    handle to handpickpeaks_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zh = zoom(handles.figure1);
set(zh, 'Enable','off');
if ~isempty(handles.data.peakpts)
  delete(handles.data.peakpts)
end
axes(handles.axes1);
[X,Y] = ginput;
X = sort(X);
Y = interp1(handles.data.samples, handles.data.raw, X, 'nearest');
% axes(handles.axes1);
hold(handles.axes1,'on');
handles.data.peakpts = plot(handles.axes1,X,Y,'.r');
set(zh,'Motion','both','Enable','on');
df = ceil(min(diff(X))*.3);
handles.data.prefilterwidth = df;
set(handles.prefilterwidth,'String',num2str(handles.data.prefilterwidth));
guidata(hObject, handles);

%--------------------------------------------------------------------------
function handles = autofindpeaks(hObject, handles)
handles.data.prefilterwidth = str2double(get(handles.prefilterwidth,'String'));
[ fx, fy ] = fitfringe(handles.data.raw, handles.data.prefilterwidth, -1, 0 );
handles.data.peakx = fx;
handles.data.peaky = fy;
F0 = 1;
if ~isempty(fx)
  V = polyfit(handles.data.peakx*1e-3,handles.data.peaky,3);
  handles.data.Y = [ handles.data.X V F0 ];
end
guidata(hObject, handles);
update_Y_to_fig(handles);

%--------------------------------------------------------------------------
% --- Executes on button press in findpeaks_btn.
function findpeaks_btn_Callback(hObject, ~, handles)
% hObject    handle to findpeaks_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = autofindpeaks(hObject, handles);
execute_level(hObject, handles);

%--------------------------------------------------------------------------
% --- Executes on button press in next_btn.
function next_btn_Callback(hObject, ~, handles)
% hObject    handle to next_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.data.level
  case 1
    handles.data.level = 2;
    guidata(hObject, handles);
    handles = setup_level(hObject, handles);
    execute_level(hObject, handles);
  case 2
    handles.data.level = 3;
    handles.data.fitpass = 1;
    handles.data.Y(1:7) = handles.data.X;
    guidata(hObject, handles);
    handles = setup_level(hObject, handles);
    execute_level(hObject, handles);
  case 3
    % Write out current fit
    guidata(hObject, handles);
    if next_scannum_file(hObject, handles)
      set(handles.next_btn,'enable','off');
      set(handles.reiterate_btn,'enable','off');
      set(handles.back_btn,'enable','off');
    end
end

% --- Executes on button press in back_btn.
function back_btn_Callback(hObject, ~, handles)
% hObject    handle to back_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.data.level > 1
  handles.data.level = handles.data.level - 1;
  if handles.data.level == 2
    if isempty(handles.data.peakx)
      handles.data.level = 1;
    end
%     if ~isempty(handles.data.Y)
%       handles.data.X = handles.data.Y(1:7);
%     end
  end
  handles.data.startlevel = handles.data.level;
  guidata(hObject,handles);
  handles = setup_level(hObject,handles);
  execute_level(hObject,handles);
end

%--------------------------------------------------------------------------
% --- Executes on button press in reiterate_btn.
function reiterate_btn_Callback(hObject, ~, handles)
% hObject    handle to reiterate_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
execute_level(hObject, handles);

%--------------------------------------------------------------------------
function handles = setup_level( hObject, handles )
if handles.data.level >= handles.data.startlevel
  switch handles.data.level
    case 1
      % hide level 2 and 3 elements
      set(handles.reiterate_btn,'enable','off');
      set(handles.back_btn,'enable','off');
      set(handles.next_btn,'enable','off','String','Next');
      set(handles.peakfit_panel,'visible','on');
      set(handles.fullfit_panel,'visible','off');
      set(handles.axes2,'visible','off');
      cla(handles.axes2);
      set(handles.axes3,'visible','off');
      cla(handles.axes3);
      set(handles.Pause,'visible','off');
      set(handles.Slider,'visible','off');
      set(handles.tc_select,'visible','off');
      % expose level 1 elements
      set(handles.peakdet_panel,'visible','on');
      set([handles.X1 handles.X2 handles.X3 handles.X4 handles.X6 ...
           handles.X5 handles.X7], ...
          'BackgroundColor',handles.data.Xcolor.full);
      cla(handles.axes1);
      plot(handles.axes1, handles.data.samples, handles.data.raw);
      handles.data.peakpts = [];
      guidata(hObject,handles);
    case 2
      set(handles.axes2,'visible','on');
      cla(handles.axes2);
      set(handles.axes3,'visible','off');
      cla(handles.axes3);
      set(handles.reiterate_btn,'enable','off');
      set(handles.back_btn,'enable','on');
      set(handles.next_btn,'enable','on','String','Next');
      set(handles.peakfit_panel,'visible','on');
      set(handles.peakfit_panel,'Title','Preliminary Fit to Peaks');
      set(handles.Slider,'visible','on');
      set([handles.X1 handles.X2 handles.X3 handles.X4 handles.X6], ...
          'BackgroundColor',handles.data.Xcolor.fit);
      if get(handles.dblexp,'Value') == 0
          set(handles.tc_select,'SelectedObject',handles.tc1_enbl);
          set(handles.tc_select,'visible','off');
      else
          set(handles.tc_select,'visible','on');
      end
      set(handles.X5,'BackgroundColor',handles.data.Xcolor.full);
      set(handles.X7,'BackgroundColor',handles.data.Xcolor.full);
      set(handles.peakdet_panel,'visible','off');
      set(handles.fullfit_panel,'visible','off');
      set(handles.Pause,'visible','off');
      % update_X_to_fig(handles);
    case 3
      set(handles.Slider,'visible','off');
      set(handles.tc_select,'visible','off');
      set(handles.reiterate_btn,'enable','on');
      set(handles.back_btn,'enable','on');
      set(handles.next_btn,'String','Next');
      if handles.data.index < length(handles.data.scans)
        set(handles.next_btn,'enable','on');
      else
        set(handles.next_btn,'enable','off');
      end
      set([handles.X1 handles.X2 handles.X3 handles.X4 handles.X6 ...
           handles.X5 handles.X7], ...
          'BackgroundColor',handles.data.Xcolor.full);
      set(handles.peakfit_panel,'visible','on');
      set(handles.peakfit_panel,'Title','Full Fit');
      set(handles.peakdet_panel,'visible','off');
      set(handles.fullfit_panel,'visible','on');
      set(handles.Pause,'Value',1,'visible','on');
      update_Y_to_fig(handles);
  end
end

%--------------------------------------------------------------------------
function all_done = next_scannum_file(hObject, handles )
% handles may be modified as a side effect.
while 1
  if handles.data.index > 0 && handles.data.figerr >= 0 && ...
      handles.data.figerr < handles.data.threshold
    handles.data.Xlast = handles.data.X;
    handles.data.Ylast = handles.data.Y;
    i = handles.data.index;
    if handles.data.saveall
      fprintf( handles.data.ofd, ...
        '%d %.2f %.1f %d %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g %.7g\n', ...
        handles.data.scans(i), handles.data.P(i), ...
        handles.data.T(i), ...
        handles.data.samples(1), handles.data.Y(1:12), ...
        handles.data.fitpass, handles.data.figerr );
    else
      fprintf( handles.data.ofd, ...
        '%d %.2f %.1f %d %.7g %.7g %.7g %.7g %.7g %.7g %.7g\n', ...
        handles.data.scans(i), handles.data.P(i), ...
        handles.data.T(i), ...
        handles.data.samples(1), handles.data.Y(1:7) );
    end
  end
  if handles.data.index >= length(handles.data.scans)
    all_done = 1;
    fclose(handles.data.ofd);
    handles.data.ofd = [];
    guidata(hObject, handles);
    return;
  else
    all_done = 0;
  end
  handles.data.index = handles.data.index+1;
  handles.data.fitpass = 1;
  scannum = handles.data.scans(handles.data.index);
  fe = loadbin(mlf_path(handles.data.ScanNumdir, scannum));
  if ~isempty(fe) && size(fe,2) >= 2 && size(fe,1) >= max(handles.data.samples)
    set(handles.SCANNUM,'String',num2str(scannum));
    if isempty(handles.data.vZ)
      handles.data.raw = fe(handles.data.samples,2);
    else
      handles.data.raw = fe(handles.data.samples,2)-mean(fe(handles.data.vZ,2));
    end
    handles.data.peakx = [];
    if handles.data.startlevel ~= handles.data.level
      handles.data.level = handles.data.startlevel;
      handles = setup_level(hObject, handles);
    end
    guidata(hObject,handles);
    if execute_level(hObject,handles);
      return;
    end
    handles = guidata(hObject);
  else
    handles.data.figerr = -1;
    handles.data.figerrs(handles.data.index) = nan;
    % handles.data.passes(handles.data.index) = 4; % total failure
    handles.data.passes = update_passes(handles,4);
    guidata(hObject,handles);
    fprintf(1, 'Error reading scannum file %d\n', scannum );
  end
end

%--------------------------------------------------------------------------
function passes = update_passes( handles, npasses )
passes = handles.data.passes;
set(handles.axes3,'units','pixels');
ax3pos = get(handles.axes3,'position');
ppp = ceil(length(passes)/ax3pos(3));
if ppp > 1
  if handles.data.index < ppp
    i = 1:handles.data.index;
  else
    i = [1-ppp:0]+handles.data.index;
  end
  if npasses > max(passes(i))
    passes(i) = npasses;
  else
    passes(handles.data.index) = npasses;
  end
else
  passes(handles.data.index) = npasses;
end

%--------------------------------------------------------------------------
function interact = execute_level(hObject, handles)
% handles has been initialized, and now we do what we need to do at the
% current operating level. This function should return a non-zero value
% if the GUI should interact. It should return zero only if the fit for
% the current scannum number is complete and we want to advance to the next
% scan.
% handles may be modified as a side effect.
interact = 1;
while 1
  switch handles.data.level
    case 1
      if handles.data.startlevel <= 1
        % axes(handles.axes1);
        ttl = '';
        if ~isempty(handles.data.peakx)
          if ~isempty(handles.data.peakpts)
            delete(handles.data.peakpts);
          end
          hold(handles.axes1, 'on');
          handles.data.peakpts = ...
            plot(handles.axes1, ...
                handles.data.peakx+handles.data.samples(1)-1, ...
                handles.data.peaky, '*m');
          hold(handles.axes1, 'off');
          set(handles.next_btn,'enable','on');
          dfx = diff(handles.data.peakx);
          dsfr = dfx(2:end)./dfx(1:end-1);
          ttl = sprintf('dsfr range: %.2f to %.2f', min(dsfr), max(dsfr));
          if any(dsfr > 1.8 | dsfr < .5)
            ttl = [ ttl ': Exceeds normal range' ];
          end
        else
          set(handles.next_btn,'enable','off');
        end
        xlabel(handles.axes1,'Sample');
        title(handles.axes1,ttl);
        zh = zoom(handles.figure1);
        set(zh,'Enable','on');
        guidata(hObject, handles);
        return;
      else
        handles = autofindpeaks(hObject, handles);
        handles.data.level = 2;
      end
    case 2
      %### Should have X displayed at level 1, I guess
      % handles = update_X_from_fig(handles);
      if isempty(handles.data.peakx)
        % We found no peaks at all, so we want to set up to skip
        % at level 3
        handles.data.level = 3;
        handles.data.figerr = -1;
        handles.data.figerrs(handles.data.index) = nan;
        handles.data.fitpass = 4;
      else
        fn = (1:length(handles.data.peakx))';
        px = handles.data.peakx*1e-3;
        XX = [ ones(size(px)) px px.*px ];
        set_fitting(handles, 2);
        dblexp = get(handles.dblexp,'Value');
        if dblexp
            M = [XX exp(-px/handles.data.X(5)) exp(-px/handles.data.X(7))];
            ffX = M\fn;
            handles.data.X(1:4) = ffX(1:4);
            handles.data.X(6) = ffX(5);
            ffn = etln_evalJ(handles.data.X, handles.data.rxs);
        else
            M = [XX exp(-px/handles.data.X(5))];
            ffX = M\fn;
            handles.data.X(1:4) = ffX(1:4);
            handles.data.X(6) = 0;
            ffn = etln_evalJ(handles.data.X(1:5), handles.data.rxs);
        end
        update_X_to_fig(handles);
        guidata(hObject,handles);
        if handles.data.startlevel <= 2
          fnm = etln_evalJ(handles.data.X, handles.data.peakx*1e-3);
          % axes(handles.axes1);
          cla(handles.axes1);
          handles.data.peakpts = [];
          x = handles.data.peakx+handles.data.samples(1)-1;
          plot(handles.axes2, x, (fn-fnm)*100, '*' );
          title(handles.axes2, 'Residual as percent of a fringe');
          % xlabel(handles.axes2, 'Sample');
          plot(handles.axes1, x, fn, '*', handles.data.samples, ffn);
          % xlabel(handles.axes1, 'Sample');
          ylabel(handles.axes1, 'Fringe Number');
          zh = zoom(handles.figure1);
          set(zh,'Enable','on');
          return;
        else
          handles.data.Y(1:7) = handles.data.X;
          handles.data.level = 3;
        end
      end
    case 3
      handles.data.startlevel = 3;
      % handles = update_Y_from_fig(handles);
      handles.data.threshold = str2num(get(handles.threshold,'String'));
      etln = handles.data.raw;
      % [ Y, resnorm,residual,exitflag,output ] = ...
      set_fitting(handles,1);
      if max(etln) - min(etln) > .1
        dblexp = get(handles.dblexp,'Value');
        if dblexp
          fY = handles.data.Y;
        else
          fY = handles.data.Y([1:5 8:12]);
        end
        fY = lsqcurvefit(@etln_evalJ, fY, ...
          handles.data.rxs, etln, ...
          [], [], handles.data.Op);
        if dblexp
          handles.data.Y = fY;
        else
          handles.data.Y([1:5 8:12]) = fY;
          handles.data.Y(6) = 0;
        end
        set_fitting(handles,0);
        % axes(handles.axes2)
        if ~any(isnan(handles.data.Y))
          emdl = etln_evalJ(handles.data.Y,handles.data.rxs);
          P = polyval(handles.data.Y(8:11),handles.data.rxs);
          handles.data.figerr = std(etln-emdl)/(max(etln)-min(etln));
          handles.data.figerrs(handles.data.index) = handles.data.figerr;
          cla(handles.axes2);
          plot(handles.axes2, handles.data.indexes,handles.data.figerrs);
          set(handles.axes2,'XTickLabel',[],'YTickLabel',[], ...
            'xlim', [1 length(handles.data.scans)+1], ...
            'ylim', [0 handles.data.threshold], 'visible','on');
          title(handles.axes2, sprintf('Relative Error: %.2g', ...
              handles.data.figerr));
        else
          handles.data.figerr = -1;
          handles.data.figerrs(handles.data.index) = nan;
          title(handles.axes2, 'Fit failed');
        end
      else
        handles.data.figerr = -1;
        handles.data.figerrs(handles.data.index) = nan;
        handles.data.fitpass = 4;
      end
      handles.data.passes = update_passes(handles,handles.data.fitpass);
      % handles.data.passes(handles.data.index) = handles.data.fitpass;
      guidata(hObject,handles);
      update_Y_to_fig(handles);
      % axes(handles.axes1);
      cla(handles.axes1);
      handles.data.peakpts = [];
      if handles.data.figerr >= 0
        plot(handles.axes1, handles.data.samples, etln, 'g', ...
          handles.data.samples, emdl, 'b', ...
          handles.data.samples, P, 'r');
      else
        plot(handles.axes1, handles.data.samples, etln, 'g' );
      end
      xlabel(handles.axes1, 'Sample');

      % axes(handles.axes3);
      cla(handles.axes3);
      colors = [ 1 1 1;
        0 1 0;
        .5 1 0;
        1 1 0;
        1 .5 0;
        1 0 0 ];
      C = zeros(1, length(handles.data.passes), 3);
      % C(1,:,:) = colors(handles.data.passes+1,:);
      C(1,:,:) = interp1([0; 1; 4.5], [1 1 1; 0 1 0; 1 0 0], handles.data.passes);
      image(C,'Parent',handles.axes3);
      %plot(handles.data.indexes,handles.data.passes);
      set(handles.axes3,'XTickLabel',[],'YTickLabel',[],'Visible','on',...
        'XTick',[],'YTick',[]);
      dopause = get(handles.Pause, 'Value');
      % dopause values are: 1: stop always, 2: stop on failure, 3: skip on
      % failure
      if handles.data.figerr >=0 && ...
              handles.data.figerr < handles.data.threshold && ...
              handles.data.Y(12) > 0
        set(handles.next_btn,'String','Next');
        interact = dopause == 1;
        return;
      elseif handles.data.Y(12) < 0 && ...
           handles.data.figerr >=0 && ...
           handles.data.figerr < handles.data.threshold && ...
           handles.data.fitpass == floor(handles.data.fitpass)
        handles.data.Y(12) = -handles.data.Y(12);
        handles.data.Y(8:11) = handles.data.Y(8:11)/(1 - handles.data.Y(12));
        handles.data.Y(1) = handles.data.Y(1) + .5;
        handles.data.fitpass = handles.data.fitpass + .5;
      elseif handles.data.figerr >= 0 && ...
              handles.data.Y(12) < 0 && ...
              handles.data.fitpass == floor(handles.data.fitpass)
        handles.data.Y(12) = -handles.data.Y(12);
        handles.data.fitpass = handles.data.fitpass + .5;
      elseif handles.data.fitpass < 2 && isempty(handles.data.peakx) && ...
          ~isempty(handles.data.Ylast) && ...
          ~any(isnan(handles.data.Ylast(1:7)))
        handles.data.X = handles.data.Ylast(1:7);
        handles.data.fitpass = 2;
      elseif handles.data.fitpass < 3 && ~isempty(handles.data.Xlast);
        handles.data.X = handles.data.Xlast;
        handles.data.fitpass = 3;
      elseif handles.data.fitpass < 4
        handles.data.X = handles.data.Xdflt;
        handles.data.fitpass = 4;
      else
        set(handles.next_btn,'String','Skip');
        if dopause > 2
          interact = 0;
        else
          interact = 1;
        end
        return;
      end
      if handles.data.fitpass ~= floor(handles.data.fitpass)
          handles.data.level = 1;
      end
      guidata(hObject,handles);
      update_X_to_fig(handles);
  end
end

%--------------------------------------------------------------------------
function update_X_to_fig(handles)
set(handles.X1, 'String', num2str(handles.data.X(1)));
set(handles.X2, 'String', num2str(handles.data.X(2)));
set(handles.X3, 'String', num2str(handles.data.X(3)));
set(handles.X4, 'String', num2str(handles.data.X(4)));
set(handles.X5, 'String', num2str(handles.data.X(5)));
set(handles.X6, 'String', num2str(handles.data.X(6)));
set(handles.X7, 'String', num2str(handles.data.X(7)));

%--------------------------------------------------------------------------
function update_Y_to_fig(handles)
handles.data.X = handles.data.Y(1:7);
update_X_to_fig(handles);
set(handles.X8, 'String', num2str(handles.data.Y(8)));
set(handles.X9, 'String', num2str(handles.data.Y(9)));
set(handles.X10, 'String', num2str(handles.data.Y(10)));
set(handles.X11, 'String', num2str(handles.data.Y(11)));
set(handles.X12, 'String', num2str(handles.data.Y(12)));

% function handles = update_X_from_fig(handles)
% X(1) = str2double(get(handles.X1,'String'));
% X(2) = str2double(get(handles.X2,'String'));
% X(3) = str2double(get(handles.X3,'String'));
% X(4) = str2double(get(handles.X4,'String'));
% X(5) = str2double(get(handles.X5,'String'));
% X(6) = str2double(get(handles.X6,'String'));
% X(7) = str2double(get(handles.X7,'String'));
% handles.data.X = X;

% function handles = update_Y_from_fig(handles)
% handles = update_X_from_fig(handles);
% Y = handles.data.X;
% Y(8) = str2double(get(handles.X8,'String'));
% Y(9) = str2double(get(handles.X9,'String'));
% Y(10) = str2double(get(handles.X10,'String'));
% Y(11) = str2double(get(handles.X11,'String'));
% Y(12) = str2double(get(handles.X12,'String'));
% handles.data.Y = Y;

% --- Executes on button press in Pause.
function Pause_Callback(~, ~, ~)
% hObject    handle to Pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Pause


function set_fitting(handles, turn_on )
if turn_on == 1
  set(handles.Fitting,'visible','on');
  set(handles.next_btn,'enable','off');
  set(handles.reiterate_btn,'enable','off');
  set(handles.back_btn,'enable','off');
  drawnow;
elseif turn_on == 0
  set(handles.Fitting,'visible','off');
  set(handles.next_btn,'enable','on');
  set(handles.reiterate_btn,'enable','on');
  set(handles.back_btn,'enable','on');
elseif turn_on == 2
  set(handles.Fitting,'visible','off');
  set(handles.next_btn,'enable','on');
  set(handles.reiterate_btn,'enable','off');
  set(handles.back_btn,'enable','on');
end

function close_request_fcn(~, ~, handles)
if ~isempty(handles.data.ofd)
  fclose(handles.data.ofd);
end
warning('on','MATLAB:rankDeficientMatrix');
delete(handles.figure1);

% --- Executes on button press in defaults_btn.
function defaults_btn_Callback(hObject, ~, handles)
if handles.data.level == 3
    handles.data.X = handles.data.Y(1:7);
    guidata(hObject,handles);
end
waveform = handles.data.wv.Name;
save_waveform_params( waveform, 'threshold', handles.data.threshold, ...
  'prefilterwidth', handles.data.prefilterwidth,'X', handles.data.X );

% --- Executes on button press in dblexp.
function dblexp_Callback(hObject, ~, handles)
dblexp = get(hObject,'Value');
if dblexp
  set(handles.dblexp1,'Visible','on');
  set(handles.dblexp2,'Visible','on');
  set(handles.dblexp3,'Visible','on');
  set(handles.X6, 'Visible', 'on');
  set(handles.X7, 'Visible', 'on');
  set(handles.tc_select,'Visible','on');
  set(handles.tc_select,'SelectedObject',handles.tc2_enbl);
else
  set(handles.dblexp1,'Visible','off');
  set(handles.dblexp2,'Visible','off');
  set(handles.dblexp3,'Visible','off');
  set(handles.X6, 'Visible', 'off');
  set(handles.X7, 'Visible', 'off');
  set(handles.tc_select,'Visible','off');
  set(handles.tc_select,'SelectedObject',handles.tc1_enbl);
end
if handles.data.level > 1
    execute_level(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function peakfit_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peakfit_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function newval = apply_step(oldval, step, handle)
newval = oldval + step;
if newval < 0
    newval = oldval;
else
    set(handle,'String',num2str(newval));
end

% --- There is probably a slicker way to do this.
function Slider_Motion(hObject, ~, handles)
% fprintf(1,'Slider_Listener event\n');
if handles.data.polling
    return;
end
handles.data.polling = 1;
guidata(hObject,handles);
while 1
    handles = guidata(hObject);
    if ~handles.data.polling
        return;
    end
    Val = get(handles.Slider,'Value');
    stepsize = 10^(abs(Val)-3);
    step = sign(Val)*stepsize;
    Movee = get(handles.tc_select,'SelectedObject'); % handles.data.SliderSelected;
    if Movee == handles.tc1_enbl
        handles.data.X(5) = ...
            apply_step(handles.data.X(5), step, handles.X5);
    elseif Movee == handles.tc2_enbl
        handles.data.X(7) = ...
            apply_step(handles.data.X(7), step, handles.X7);
    end
    guidata(hObject,handles);
    execute_level(hObject, handles);
    pause(.1);
end


% --- Executes on slider movement.
function Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.Slider,'Value',0);
handles.data.polling = 0;
guidata(hObject,handles);


% % --- If Enable == 'on', executes on mouse press in 5 pixel border.
% % --- Otherwise, executes on mouse press in 5 pixel border or over X5.
% function X57_ButtonDownFcn(hObject, eventdata, handles)
% % hObject    handle to X5 or X7 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% handles.data.SliderSelected = hObject;
% guidata(hObject, handles);
% handles = setup_level(hObject, handles);
% execute_level(hObject, handles);
