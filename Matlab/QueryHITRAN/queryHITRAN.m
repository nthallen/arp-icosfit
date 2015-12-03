function queryHITRAN( varargin )

%queryHITRAN; start up the gui
%You should edit this file and replace user and pword with your MySQL
%username and password.
%In order to use this program the HITRAN database needs to be in MyQSL
%database. 
%From the command line
%$ gzip -d HITRAN.sql.gz
%$ mysql --username 'sayres' --password < HITRAN.sql
%
%Depending on your platform you may need to tell matlab where the mysql
%driver is. For example on linux: javaaddpath /usr/share/java/mysql-connector-java.jar
%add to file. Search matlab help on mysql and adding JDBC driver. 

DB='HITRAN';
HITRAN_Table='HITRAN2012';
MOLEC_INFO_Table='HITRAN_Molecules';
user='sayres'; %change this to your MyQSL user name
pword='mysqlpw'; %change this to your MyQSL passowrd

if nargin == 0
line_obj.DB=DB;
line_obj.HITRAN_Table=HITRAN_Table;
line_obj.MOLEC_INFO_Table=MOLEC_INFO_Table;
line_obj.user=user;
line_obj.pword=pword;

if isempty(line_obj.user)
    mysql('open');
else
    conn=database('HITRAN',user,pword,'Vendor','MySQL');
end
setdbprefs('DataReturnFormat','cellarray')
[fields]=fetch(conn,['SHOW FIELDS FROM ' HITRAN_Table]);
fields=fields(:,1);
[line_obj.tables]=fetch(conn,'SHOW TABLES');
setdbprefs('DataReturnFormat','structure')
[line_obj.HITRAN]=fetch(conn,['select molec,iso,isotopologue,weight,frac_abun,name,texname from ' MOLEC_INFO_Table]);
close(conn);
line_obj.fields=fields;
line_obj.sel_fig = 0;
line_obj.query_results = [];
for i=1:length(line_obj.HITRAN.molec)
    eval(['line_obj.select.m' num2str(line_obj.HITRAN.molec(i)) '.i' num2str(line_obj.HITRAN.iso(i)) '=1;']);
end
for i=1:max(line_obj.HITRAN.molec)
    eval(['line_obj.select.m' num2str(i) '.i0 = 0;']);
end
%[trop, strat, user]
line_obj.mixing_ratios=[...
    10000e-6, 5e-6, 100e-6;... %H2O
    388e-6, 388e-6, 388e-6;... %CO2
    20e-9, 1e-6, 80e-9;... %O3
    0.3e-6, 0.3e-6, 0.3e-6;... %N2O
    20e-9, 20e-9, 20e-9;... %CO
    1.7e-6, 1.5e-6, 1.7e-6;... %CH4
    .2, .2, .2;... %O2
    1e-12,1e-9,0;... %NO
    0,0,0;... %SO2
    1e-12,5e-9,0;... %NO2
    0,0,0;... %NH3
    5e-12,1e-8,1e-10,;... %HNO3
    1e-13,2e-12,0,;... %OH
    0,0,0,;... %HF
    0,8e-11,0;... %HCl
    0,0,0;... %HBr
    0,0,0;... %HI
    0,2e-10,0;... %ClO
    0,0,0;... %OCS
    0,0,0;... %H2CO
    0,0,0;... %HOCl
    .78,.78,.78;... %N2
    0,0,0;... %HCN
    0,0,0;... %CH3Cl
    1e-8,1e-8,0;... %H2O2
    11e-9,0,0;... %C2H2
    11e-9,0,0;... %C2H6
    0,0,0;... %PH3
    0,5e-11,0;... %COF2
    0,0,0;... %SF6
    0,0,0;... %H2S
    0,0,0;... %HCOOH
    1e-11,1e-11,0;... %HO2
    0,0,0;... %O
    0,5e-11,0;... %ClONO2
    0,0,0;... %NO+
    0,0,0;... %HOBr
    11e-9,0,0;... %C2H4
    0,0,0;... %CH3OH
    0,0,0;... %CH3Br
    0,0,0;... %CH3CN
    0,0,0;... %CF4
    0,0,0;... %C4H2
    0,0,0;... %HC3N
    0,0,0;... %H2
    0,0,0;... %CS
    0,0,0;... %SO3
    ];
cmap=colormap('lines');
cmap(8:15,:)=[255 153 0;
    153 153 255;
    153 153 0;
    246 137 249;
    255 102 102;
    110 187 154;
    187 42 236;
    230 37 162;]/255;
line_obj.cmap=cmap;
%set up output field parameters
f = figure('visible','off','Name','Query HITRAN Database Program','Numbertitle','off');
set(f,'UserData',line_obj,'tag','queryHITRAN')
hp(1)=uipanel('Title','Select Parameters','FontSize',12,...
                'Position',[.02 .1 .6 .85],'Unit','normalized');
h = zeros(length(fields),1);
en = zeros(length(fields),1);
for i = 1:length(fields)
    if strcmp(fields(i,:),'molec') || ...
            strcmp(fields(i,:),'iso') || ...
            strcmp(fields(i,:),'wavenumber') || ...
            strcmp(fields(i,:),'intensity') || ...
            strcmp(fields(i,:),'gamma_air') || ...
            strcmp(fields(i,:),'n_air') || ...
            strcmp(fields(i,:),'delta_air') || ...
            strcmp(fields(i,:),'Ierr') || ...
            strcmp(fields(i,:),'LS_energy')
        en(i)=1;
    end
end
height = 0;
width = 0;
for i=length(fields):-1:1
    name = fields(i,:);
    h(i) = uicontrol( hp(1), 'Style', 'checkbox','tag','parameters','UserData',i,'string',name,'value',en(i));
    pos = get(h(i),'extent');
    height = height + pos(4);
    width = max(width,pos(3));
    set(h(i),'position',[ 220 height pos(3)+20 pos(4)]);
end
height = height+20;
uicontrol(hp(1),'style','text','string','Output Parameters','position',[220 height 130 15]);
maxheight=height+20;
height=height-20;
%set up select parameters
hs(1)=uicontrol(hp(1),'Style','text','string','Wavenumber Range','position',[10 height 150 20]);
pos=get(hs(1),'extent'); height=height-pos(4);
hs(2)=uicontrol(hp(1),'Style','text','string','min:','position',[10 height 30 20]);
hs(3)=uicontrol(hp(1),'Style','edit','tag','min_wn','string','','position',[ 40 height 50 20]);
hs(4)=uicontrol(hp(1),'Style','text','tag','wn','string','cm-1','position',[90 height 50 20]);
pos=get(hs(2),'extent'); height=height-pos(4);
hs(5)=uicontrol(hp(1),'Style','text','string','max:','position',[10 height 30 20]);
hs(6)=uicontrol(hp(1),'Style','edit','tag','max_wn','string','','position',[ 40 height 50 20]);
hs(7)=uicontrol(hp(1),'Style','text','tag','wn','string','cm-1','position',[90 height 50 20]);
uicontrol(hp(1),'Style','pushbutton','string',sprintf('Switch\n cm-1/um'),'Callback','queryHITRAN(''switch_wn'')','position',[140 height 60 40]);
pos=get(hs(7),'extent'); height=height-pos(4)-20;
hs(8)=uicontrol(hp(1),'Style','pushbutton','string','Select Molecules','tag','sel_molec','Callback','queryHITRAN(''select_molecule'')','position',[10 height 120 25]);
pos=get(hs(7),'extent'); height=height-pos(4)-10;
hs(9)=uicontrol(hp(1),'Style','text','string','Cutoff:','position',[10 height 50 20]);
hs(10)=uicontrol(hp(1),'Style','edit','tag','cutoff','string','','position',[ 60 height 60 20]);
pos=get(hs(9),'extent'); height=height-pos(4)-100;
uicontrol(hp(1),'Style','pushbutton','string','Run Query','Callback','queryHITRAN(''run_query'')','position',[30 height 100 60]);
hs(11)=uicontrol(hp(1),'Style','text','string',sprintf('Database = %s',line_obj.DB));
pos=get(hs(11),'extent'); height=height-pos(4)-50;
set(hs(11),'Position',[10 height pos(3) pos(4)]);
hs(12)=uicontrol(hp(1),'Style','text','string','HITRAN Table = ');
pos=get(hs(12),'extent'); height=height-pos(4);
set(hs(12),'Position',[10 height pos(3) pos(4)]);
itable=find(strcmp(line_obj.tables,line_obj.HITRAN_Table));
hs(13)=uicontrol(hp(1),'Style','listbox','string',[line_obj.tables],'max',1,'min',1,'value',itable,'Callback','queryHITRAN(''change_table'')','tag','table_name');
xstart=10+pos(3);
pos=get(hs(13),'extent');
set(hs(13),'Position',[xstart height pos(3)+20 pos(4)]);



%set up select save options
hp(2)=uipanel('Title','Save Functions','FontSize',12,...
                'Position',[.65 .71 .33 .27]);
height=90;
hf(1)=uicontrol(hp(2),'Style','text','string','Fitline File Name:','position',[10 height 115 20]);
pos=get(hf(1),'extent'); height=height-pos(4);
hf(2)=uicontrol(hp(2),'Style','edit','tag','filename','string','fitline.dat','position',[10 height 150 20]);
pos=get(hf(1),'extent'); height=height-pos(4)-10;
hf(3)=uicontrol(hp(2),'Style','pushbutton','string','Create Fitline File','Callback','queryHITRAN(''save_queryHITRAN'')','position',[ 10 height 150 25]);
pos=get(hf(2),'extent'); height=height-pos(4)-10;
hf(5)=uicontrol(hp(2),'Style','pushbutton','string','Save to Workspace','Callback','queryHITRAN(''save_workspace'')','position',[10 height 150 25]);

%set up select plot options
hp(3)=uipanel('Title','Plot Functions','FontSize',12,...
                'Position',[.65 .02 .33 .66]);
set(hp(3),'Units','pixels');
pos=get(hp(3),'Position');
height=pos(4)-10;
set(hp(3),'Units','normalized');
hf(3)=uicontrol(hp(3),'Style','pushbutton','string','Intensity Plot','Callback','queryHITRAN(''plot_intensity'')','position',[ 10 height 150 25]);
pos=get(hf(2),'extent'); height=height-pos(4)-40;

hf(20)=uibuttongroup('Parent',hp(3),'Position',[0 .8 1 .1],'tag','intensity_method');
ub(1)=uicontrol(hf(20),'Style','radio','string','Intensity','Position',[5 2 70 20]);
ub(2)=uicontrol(hf(20),'Style','radio','string','Intensity*[Conc.]','Position',[80 2 120 20]);

hf(5)=uicontrol(hp(3),'Style','pushbutton','string','Absorption Spectra','Callback','queryHITRAN(''plot_transmission'')','position',[10 height 150 25]);

hf(6)=uibuttongroup('Parent',hp(3),'Position',[0 .60 1 .1],'tag','absorb_method');
ub(1)=uicontrol(hf(6),'Style','radio','string','Direct','Position',[5 2 60 20]);
ub(2)=uicontrol(hf(6),'Style','radio','string','ICOS','Position',[70 2 60 20]);

pos=get(hf(5),'extent'); height=height-pos(4)-35;
hf(7)=uicontrol(hp(3),'Style','text','string','Cell Properties');
pos=get(hf(7),'extent');
set(hf(7),'Position',[40 height pos(3) pos(4)]);
height=height-pos(4);
hf(8)=uicontrol(hp(3),'Style','text','string','Length:');
pos=get(hf(8),'extent');
set(hf(8),'Position',[5 height pos(3) pos(4)]);
hf(9)=uicontrol(hp(3),'Style','edit','tag','length','string','','Position',[5+pos(3) height 40 pos(4)]);
hf(10)=uicontrol(hp(3),'Style','text','string','cm');
pos2=get(hf(10),'extent');
set(hf(10),'Position',[5+pos(3)+40 height pos2(3) pos2(4)]);
xstart=5+pos(3)+40+pos2(3)+5;
hf(8)=uicontrol(hp(3),'Style','text','string','R:');
pos=get(hf(8),'extent');
set(hf(8),'Position',[xstart height pos(3) pos(4)]);
hf(9)=uicontrol(hp(3),'Style','edit','tag','R','string','','Position',[xstart+pos(3) height 40 pos(4)]);
hf(10)=uicontrol(hp(3),'Style','text','string','ppm');
pos2=get(hf(10),'extent');
set(hf(10),'Position',[xstart+pos(3)+40 height pos2(3) pos2(4)]);
height=height-pos(4)-5;
hf(11)=uicontrol(hp(3),'Style','text','string','Temp:');
pos=get(hf(11),'extent');
set(hf(11),'Position',[5 height pos(3) pos(4)]);
hf(12)=uicontrol(hp(3),'Style','edit','tag','temp','string','','Position',[5+pos(3) height 40 pos(4)]);
hf(13)=uicontrol(hp(3),'Style','text','string','K');
pos2=get(hf(13),'extent');
set(hf(13),'Position',[5+pos(3)+40 height pos2(3) pos2(4)]);
xstart=5+pos(3)+40+pos2(3)+5;
hf(14)=uicontrol(hp(3),'Style','text','string','Pres:');
pos=get(hf(14),'extent');
set(hf(14),'Position',[xstart height pos(3) pos(4)]);
hf(15)=uicontrol(hp(3),'Style','edit','tag','pres','string','','Position',[xstart+pos(3) height 40 pos(4)]);
hf(16)=uicontrol(hp(3),'Style','text','string','Torr');
pos2=get(hf(16),'extent');
set(hf(16),'Position',[xstart+pos(3)+40 height pos2(3) pos2(4)]);

hf(17)=uibuttongroup('Parent',hp(3),'Position',[0 .17 1 .15],'tag','profile');
uicontrol(hf(17),'Style','text','string','Mixing Ratio Profile','Position',[2 20 200 20]);
ub(1)=uicontrol(hf(17),'Style','radio','string','Trop','tag','1');
pos1=get(ub(1),'extent');
set(ub(1),'Position',[2 1 pos1(3)+20 pos1(4)])
ub(2)=uicontrol(hf(17),'Style','radio','string','Strat','tag','2');
pos2=get(ub(2),'extent');
set(ub(2),'Position',[2+pos1(3)+20 1 pos2(3)+20 pos2(4)])
ub(3)=uicontrol(hf(17),'Style','radio','string','User','tag','3');
pos3=get(ub(3),'extent');
set(ub(3),'Position',[2+pos1(3)+pos2(3)+40 1 pos3(3)+20 pos3(4)])
hf(18)=uicontrol(hf(17),'Style','pushbutton','string','Define','Callback','queryHITRAN(''MR_define'')','position',[2+pos1(3)+pos2(3)+60+pos3(3) 2 50 20]);

hf(19)=uibuttongroup('Parent',hp(3),'Position',[0 .005 1 .16],'tag','line_method');
ub(1)=uicontrol(hf(19),'Style','radio','string','Convolute Lines');
pos1=get(ub(1),'extent');
set(ub(1),'Position',[2 23 pos1(3)+20 pos1(4)])
ub(2)=uicontrol(hf(19),'Style','radio','string','Overlay Lines');
pos2=get(ub(2),'extent');
set(ub(2),'Position',[2 2 pos2(3)+20 pos2(4)])

te_x = 220+width;
br_x = te_x + 130;
br_xx = br_x + width;
pos = get(f,'position');
   pos(3) = br_xx + 170;
   pos(4) = maxheight + 90;
  set(f,'position', pos, 'resize', 'off', 'Menubar','none',...
    'UserData',line_obj,'visible','on');


elseif strcmp(varargin{1},'select_molecule')
  f=get(gcbo,'parent');
  f=get(f,'parent');
  line_obj = update_line_obj(f);
  if length(line_obj)
    if line_obj.sel_fig && ishandle(line_obj.sel_fig)
      set(line_obj.sel_fig,'visible','on');
      figure(line_obj.sel_fig);
    else

  ro.parent = f;
  line_obj.sel_fig = figure('visible','off','MenuBar','none',...
        'DeleteFcn','queryHITRAN(''sel_mol_deletefig'')','UserData',ro);
  set(f,'UserData',line_obj,'DeleteFcn','queryHITRAN(''deletefig'')');
  
uicontrol( line_obj.sel_fig, 'Style', 'checkbox','tag','All','string','Select All','value',0,...
    'CallBack','queryHITRAN(''select_all'')',...
    'Position',[10 350 100 20]);
uicontrol( line_obj.sel_fig, 'Style', 'checkbox','tag','Seven','string','Select First 7','value',0,...
    'CallBack','queryHITRAN(''select_seven'')',...
    'Position',[100 350 150 20]);
uicontrol( line_obj.sel_fig, 'Style', 'checkbox','tag','None','string','Select None','value',0,...
    'CallBack','queryHITRAN(''select_none'')',...
    'Position',[230 350 100 20]);
uicontrol( line_obj.sel_fig, 'Style', 'pushbutton','tag','Save','string','Save Selections',...
    'CallBack','queryHITRAN(''select_save'')',...
    'Position',[360 350 120 30]);
%Read and Create Molecule Table
  
height=320;
xstart=10;
  for i = 1:max(line_obj.HITRAN.molec)
      hc=uicontrol( line_obj.sel_fig, 'Style', 'checkbox','tag',num2str(i),'String','','Position',[xstart,height,20,20]);
      set(hc,'Value',eval(['line_obj.select.m' num2str(i) '.i0']));
      hb=uicontrol( line_obj.sel_fig, 'Style', 'pushbutton','tag',num2str(i),...
          'String',line_obj.HITRAN.name(line_obj.HITRAN.molec == i & line_obj.HITRAN.iso == 1),'CallBack','queryHITRAN(''select_iso'')');
      pos=get(hb,'extent');
      set(hb,'Position',[xstart+20,height,pos(3)+20,pos(4)]);
      
      if mod(i,6) == 0
          height=height-30;
          xstart=10;
      else
          xstart=xstart+85;
      end
  end
  set(f,'UserData',line_obj);
set(line_obj.sel_fig,'visible','on','name','Select Molecules/Isotopes','NumberTitle','off');
      drawnow; shg;
    end
  end
  return
  
elseif strcmp(varargin{1},'deletefig')
  f = gcbo;
  line_obj = get(f,'UserData');
  if line_obj.sel_fig && ishandle(line_obj.sel_fig)
    delete(line_obj.sel_fig);
  end
  return
  
elseif strcmp(varargin{1},'sel_mol_deletefig')

  return
  
elseif strcmp(varargin{1},'select_save')
   f=get(gcbo,'parent');
   ro=get(f,'UserData');
   line_obj=get(ro.parent,'UserData');
   h=findobj(f,'Style','checkbox');
   for i=1:length(h)
       mol=str2double(get(h(i),'tag'));
       if ~isempty(mol)
            eval(['line_obj.select.m' num2str(mol) '.i0=get(h(i),''Value'');']);
       end
   end
   set(ro.parent,'UserData',line_obj);
    
  return
  
elseif strcmp(varargin{1},'select_all')
  f=get(gcbo,'parent');
  h=findobj(f,'Style','checkbox');
  set(h,'Value',1)
  set(findobj(f,'tag','None'),'Value',0)
  set(findobj(f,'tag','Seven'),'Value',0)
  return
    
elseif strcmp(varargin{1},'select_seven')
 f=get(gcbo,'parent');
  h=findobj(f,'Style','checkbox');
  for i=1:length(h)
      if str2double(get(h(i),'tag'))<=7
        set(h(i),'Value',1)
      else
          set(h(i),'Value',0)
      end
  end
  set(findobj(f,'tag','None'),'Value',0)
  set(findobj(f,'tag','All'),'Value',0)
  set(findobj(f,'tag','Seven'),'Value',1)
  return
    
elseif strcmp(varargin{1},'select_none')
 f=get(gcbo,'parent');  
  h=findobj(f,'Style','checkbox');
  set(h,'Value',0)
  set(findobj(f,'tag','All'),'Value',0)
  set(findobj(f,'tag','Seven'),'Value',0)
  set(findobj(f,'tag','None'),'Value',1)
  return
    
elseif strcmp(varargin{1},'select_iso')
   f=get(gcbo,'Parent');
   ro=get(f,'UserData');
   line_obj=get(ro.parent,'UserData');
   molec_num=str2double(get(gcbo,'tag'));
   iso_num=find(line_obj.HITRAN.molec==molec_num);
   fiso=figure('visible','off','MenuBar','none',...
        'DeleteFcn','queryHITRAN(''sel_iso_deletefig'')',...
        'tag',num2str(molec_num));
   pos=get(fiso,'Position');
    height=10;
    for i=length(iso_num):-1:1
       h=uicontrol('Style', 'checkbox','tag',num2str(i),...
           'String',sprintf('%s \t\t %f',cell2mat(line_obj.HITRAN.name(iso_num(i))),line_obj.HITRAN.frac_abun(iso_num(i))));
       set(h,'Value',eval(['line_obj.select.m' num2str(molec_num) '.i' num2str(i)]));
       pos2=get(h,'extent');
       set(h,'Position',[10 height pos2(3)+20 pos2(4)]);
       height=height+25;
    end
  pos=[pos(1) pos(2) 200 height];
  set(fiso,'position', pos, 'resize', 'off', 'Menubar','none',...
    'UserData',ro,'visible','on','Name','Select Iso','NumberTitle','off');
  return
  
elseif strcmp(varargin{1},'sel_iso_deletefig')
   ro=get(gcbo,'UserData');
   line_obj=get(ro.parent,'UserData');
   h=findobj(gcbo,'Style','checkbox');
   mol=get(gcbo,'tag');
   for i=1:length(h)
       iso=str2double(get(h(i),'tag'));
       if ~isempty(iso)
            eval(['line_obj.select.m' num2str(mol) '.i' num2str(iso) '=get(h(i),''Value'');']);
       end
   end
   set(ro.parent,'UserData',line_obj);
return
    
elseif strcmp(varargin{1},'switch_wn')
    h=findobj(get(gcbo,'Parent'),'tag','wn');
    if strcmp(get(h(1),'string'),'cm-1')
        set(h,'string','um')
    elseif strcmp(get(h(1),'string'),'um')
        set(h,'string','cm-1')
    end
    h(1)=findobj(get(gcbo,'Parent'),'tag','min_wn');
    h(2)=findobj(get(gcbo,'Parent'),'tag','max_wn');
    for i=1:2
        v=str2double(get(h(i),'String'));
        v=1e4/v;
        set(h(i),'String',num2str(v));
    end
    
    return
    
elseif strcmp(varargin{1},'change_table')
    f=get(get(gcbo,'parent'),'Parent');
    line_obj = get(f,'UserData');
    h=findobj(get(gcbo,'Parent'),'tag','table_name');
    t=get(h,'Value');
    line_obj.HITRAN_Table=cell2mat(line_obj.tables(t));
    set(f,'UserData',line_obj);
    return
    
elseif strcmp(varargin{1},'run_query')
  f=get(get(gcbo,'parent'),'Parent');
  line_obj = update_line_obj(f);
  sel_molec=[];
  for i=1:max(line_obj.HITRAN.molec)
      if eval(['line_obj.select.m' num2str(i) '.i0']) == 1
          for j=1:10
              try
                if eval(['line_obj.select.m' num2str(i) '.i' num2str(j)]) == 1
                    sel_molec=[sel_molec,10*i+j];
                end
              catch
              end
          end
      end
  end
  params=vertcat(line_obj.outparams{:});
  params1=sprintf('line_obj.query_results.%s,',params{:});
  params2=sprintf('%s,',params{:});
  sel_molec=sprintf('%i,',sel_molec);
  msgbox('Running Query...','Database Query','replace');
  if isempty(line_obj.user)
    conn=database(line_obj.DB);
  else
    %mysql('open','localhost',line_obj.user,line_obj.pword);
    conn=database(line_obj.DB,line_obj.user,line_obj.pword,'Vendor','MySQL');
  end
  setdbprefs('DataReturnFormat','structure')
  %mysql(['use ' line_obj.DB])
  tic;
  %eval(['[' params1(1:end-1) ']=mysql(''select ' params2(1:end-1) ' from ' line_obj.HITRAN_Table ' where CONCAT(molec,iso) IN (' sel_molec(1:end-1) ') AND wavenumber BETWEEN ' line_obj.min_wn ' AND ' line_obj.max_wn ' AND intensity >= ' line_obj.cutoff ''');']);
  eval(['line_obj.query_results=fetch(conn,''select ' params2(1:end-1) ' from ' line_obj.HITRAN_Table ' where CONCAT(molec,iso) IN (' sel_molec(1:end-1) ') AND wavenumber BETWEEN ' line_obj.min_wn ' AND ' line_obj.max_wn ' AND intensity >= ' line_obj.cutoff ''');']);
  query_time=toc;
  lines_selected=length(line_obj.query_results.molec);
  fprintf('\nQuery Took %.2f seconds.\n Number of line retrieved: %i\n',query_time,lines_selected)
  msgbox(sprintf('Running Query...\nQuery Finished\nQuery Took %.2f seconds.\n Number of line retrieved: %i\n\n',query_time,lines_selected),'Database Query','replace');
%  mysql('close')
    close(conn)
    set(f,'UserData',line_obj);
  return

elseif strcmp(varargin{1},'save_queryHITRAN')
  f=get(get(gcbo,'parent'),'Parent');
  line_obj=get(f,'UserData');
  if isempty(line_obj.query_results)
      errordlg('You must run a Query before using Plot or Save functions','No Query Found')
      return
  end
  h=findobj(get(gcbo,'Parent'),'tag','filename');
  filename=get(h,'String');
  results=struct2cell(line_obj.query_results);
  fid=fopen(filename,'w');
  for i=1:size(results{1},1);
      for j=1:size(results,1);
          if isnumeric(results{j}(i))
            fprintf(fid,'%s ',num2str(results{j}(i)));
          else
            fprintf(fid,'%s',cell2mat(results{j}(i)));
          end
      end
      fprintf(fid,'\n');
  end
  
  return
  
elseif strcmp(varargin{1},'save_workspace')
  f=get(get(gcbo,'parent'),'Parent');
  line_obj=get(f,'UserData');
  if isempty(line_obj.query_results)
      errordlg('You must run a Query before using Plot or Save functions','No Query Found')
      return
  end
  
  return
  
elseif strcmp(varargin{1},'plot_intensity')
  f=get(get(gcbo,'parent'),'Parent');
  line_obj=get(f,'UserData');
  if isempty(line_obj.query_results)
      errordlg('You must run a Query before using Plot or Save functions','No Query Found')
      return
  end
  h=get(gcbo,'Parent');
  intensity_method=get(get(findobj(h,'tag','intensity_method'),'SelectedObject'),'String');
  profile=str2double(get(get(findobj(h,'tag','profile'),'SelectedObject'),'tag'));
  figure; hold on;
  cmap=line_obj.cmap;
  k=1;
  leg=[];
  nmolec=unique(line_obj.query_results.molec);
  min_bv=min(line_obj.query_results.intensity)*.9;
  for i=nmolec'
      niso=unique(line_obj.query_results.iso(line_obj.query_results.molec==i));
      for j=niso'
        if strcmp(intensity_method,'Intensity')
            conc=1;
        else
            conc=line_obj.mixing_ratios(i,profile);   
        end
        index=find(line_obj.query_results.iso==j & line_obj.query_results.molec==i);
        h=stem(line_obj.query_results.wavenumber(index),line_obj.query_results.intensity(index)*conc,'Marker','none','color',cmap(k,:),'linewidth',2);
        set(h,'ShowBaseLine','off','BaseValue',min_bv*conc)
        set(gca,'Yscale','log')
        leg=[leg,line_obj.HITRAN.texname(find(line_obj.HITRAN.molec==i & line_obj.HITRAN.iso==j))];
        delete(get(h,'Baseline'));
        k=k+1;
      end
      k=k+3;
  end
  legend(leg);
  box on;
  if strcmp(intensity_method,'Intensity')
        ylabel('Intensity (cm^{-1}/(molecule cm^{-2}))')
  else
        ylabel('Intensity (cm^{-1})')
  end
  xlabel('Wavenumber (cm^{-1})')
  return

  elseif strcmp(varargin{1},'MR_define')
   f=get(get(get(gcbo,'parent'),'Parent'),'Parent');
   line_obj=get(f,'UserData');
   fmr=figure('visible','off','MenuBar','none',...
        'DeleteFcn','queryHITRAN(''MR_define_deletefig'')',...
        'UserData',f,'tag','MR_define',...
        'Name','Volume Mixing Ratios','NumberTitle','off');
    pos=get(fmr,'Position');
    set(fmr,'Position',[pos(1) pos(2) pos(3) 530])
    pos=get(fmr,'Position');
    height=pos(4)-25;
    xstart=0;
    for i=1:size(line_obj.mixing_ratios,1)
       if mod(i,20)==1
        uicontrol(fmr,'Style','text','String','Molecule','Position',[xstart+5 height 100 20 ]);
        uicontrol(fmr,'Style','text','String','Troposphere','Position',[xstart+110 height 100 20 ]);
        uicontrol(fmr,'Style','text','String','Stratosphere','Position',[xstart+220 height 100 20 ]);
        uicontrol(fmr,'Style','text','String','User Defined','Position',[xstart+330 height 100 20 ]);
       end
       height=height-25;
       mol_name=line_obj.HITRAN.molec==i & line_obj.HITRAN.iso==1;
       h=uicontrol(fmr,'Style','text','string',sprintf('%s',line_obj.HITRAN.name{mol_name}),'Position',[xstart+5 height 100 20]);
       h=uicontrol(fmr,'Style','text','string',sprintf('%.2e',line_obj.mixing_ratios(i,1)),'Position',[xstart+110 height 100 20]);
       h=uicontrol(fmr,'Style','text','string',sprintf('%.2e',line_obj.mixing_ratios(i,2)),'Position',[xstart+220 height 100 20]);
       h=uicontrol(fmr,'Style','edit','tag',num2str(i),'string',sprintf('%.2e',line_obj.mixing_ratios(i,3)),'Position',[xstart+330 height 100 20]);
       if mod(i,20)==0
           xstart=xstart+440;
           height=pos(4)-25;
       end
    end

    
  pos=[pos(1) pos(2) xstart+440 530];
  set(fmr,'position', pos, 'resize', 'off', 'Menubar','none',...
    'visible','on');
  return
  
elseif strcmp(varargin{1},'MR_define_deletefig')
   f=get(gcbo,'UserData');
   line_obj=get(f,'UserData');
   h=findobj(gcbo,'Style','edit');
   for i=1:length(h)
       line_obj.mixing_ratios(str2num(get(h(i),'tag')),3)=str2num(get(h(i),'string'));
   end
   set(f,'UserData',line_obj);
return
  
elseif strcmp(varargin{1},'plot_transmission')
  f=get(get(gcbo,'parent'),'Parent');
  line_obj=get(f,'UserData');
  if isempty(line_obj.query_results)
      errordlg('You must run a Query before using Plot or Save functions','No Query Found')
      return
  end
  h=get(gcbo,'Parent');
  absorb_method=get(get(findobj(h,'tag','absorb_method'),'SelectedObject'),'String');
  line_method=get(get(findobj(h,'tag','line_method'),'SelectedObject'),'String');
  profile=str2double(get(get(findobj(h,'tag','profile'),'SelectedObject'),'tag'));
  L=str2double(get(findobj(h,'tag','length'),'string'));
  R=str2double(get(findobj(h,'tag','R'),'string'));
  R=1-R*1e-6;
  T=str2double(get(findobj(h,'tag','temp'),'string'));
  P=str2double(get(findobj(h,'tag','pres'),'string'));
  molec=line_obj.query_results.molec;
  iso=line_obj.query_results.iso;
  v0=line_obj.query_results.wavenumber;
  pcoeff=line_obj.query_results.gamma_air;
  pshift=line_obj.query_results.delta_air;
  sigma=line_obj.query_results.intensity;
  E=line_obj.query_results.LS_energy;
  n=line_obj.query_results.n_air;
  
  %Now we plot by calling other functions
  if range(v0) < 10
      dv=0.001;
  elseif range(v0) < 100
      dv=0.001;
  else 
      dv=0.01;
  end
  v=[floor(min(v0)):dv:ceil(max(v0))];
  Power=ones(1,length(v))*50000;
  M=ones(length(v0),1);
  for i=unique(molec*100+iso)'
      M(molec*100+iso==i)=line_obj.HITRAN.weight(line_obj.HITRAN.molec*100+line_obj.HITRAN.iso==i);
  end
  ppm=line_obj.mixing_ratios(molec,profile);
  
  if strcmp(line_method,'Convolute Lines')
    [v,a]=absorbance(v,L,v0,pcoeff,pshift,sigma,E,n,M,T,P,ppm);
    a=a';
  elseif strcmp(line_method,'Overlay Lines')
      k=1;
      for i=unique(molec*10+iso)'
          j=find(molec*10+iso==i);
          [v,a(:,k)]=absorbance(v,L,v0(j),pcoeff(j),pshift(j),sigma(j),E(j),n(j),M(j),T,P,ppm(j));
          k=k+1;
      end
  end
  figure; hold on; box on;
  for k=1:size(a,2)
  if strcmp(absorb_method,'Direct')
    Pout=Power.*exp(-a(:,k)');
  elseif strcmp(absorb_method,'ICOS')
    delta=1e-6;  %maximum fractional error in skew
    c=2.99792458e10;
    Ts=50e-6;
    tuningrate=50;
    dt=(mean(diff(v))/tuningrate)^-1;
    % calculate number of interations
    N = c/(2*L*dt);
    M=log10(delta)/(2*N*log10(R(1)));
    dt=1/dt;
    % assumes that vector a, etc. are linear with time. (ie. a(1)=t1,
    % a(2)=t2, ...)
    tau = L./(c*(1-R+a(:,k)'));
    laserpower=Power.*Ts*c.*tau.*(1-exp(-1*dt./tau))/(2*L);
    taufactors = exp(-1*dt./tau');
    i=[0:M];   m=length(i);   n=length(taufactors);
    B=(taufactors*ones(1,m)).^(ones(n,1)*i);
    taumatrix=spdiags(B,-[0:M],length(v),length(v))';
    Pout = (laserpower*taumatrix)';
  end
  
    Pout=Pout/max(Pout)*100;
    line(v,Pout,'linewidth',1,'color',line_obj.cmap(k,:))

    if strcmp(line_method,'Convolute Lines')
      for i=1:length(molec)
        name=line_obj.HITRAN.texname(line_obj.HITRAN.molec==molec(i) & line_obj.HITRAN.iso==iso(i));
        y=min(Pout(abs(v-v0(i)) < dv*5));
        text(v0(i),y(1),name,'Rotation',270,'FontWeight','bold');
      end
    else
        mi=unique(molec*10+iso);
        name=line_obj.HITRAN.texname(line_obj.HITRAN.molec*10+line_obj.HITRAN.iso==mi(k));
      for i=find(molec*10+iso==mi(k))'
        y=min(Pout(abs(v-v0(i)) < dv*5));
        text(v0(i),y(1),name,'Rotation',270,'FontWeight','bold','color',line_obj.cmap(k,:));
      end
    end
  end
    set(gca,'XDir','reverse')
    grid
    xlabel('Wavenumber (cm^{-1})')
    ylabel('% Transmission')
    addzoom
    orient landscape
  return

end

function line_obj = update_line_obj(f)
line_obj = get(f,'UserData');
fields = { 'min_wn','max_wn','cutoff' };
for j = 1:length(fields)
  tag = fields{j};
  h = findobj(f,'tag',tag);
  for i = 1:length(h)
    line_obj.(tag) = get(h(i),'String');
  end
end
h=findobj(f,'tag','wn');
if strcmp(get(h,'String'),'um')
   line_obj.min_wn=num2str(1e4/str2double(line_obj.min_wn));
   line_obj.max_wn=num2str(1e4/str2double(line_obj.max_wn));
end
if str2double(line_obj.min_wn)>str2double(line_obj.max_wn)
    temp=line_obj.min_wn;
    line_obj.min_wn=line_obj.max_wn;
    line_obj.max_wn=temp;
end
h=findobj(f,'tag','parameters','value',1);
line_obj.outparams=get(h,'String');

set(f, 'UserData', line_obj );
return;

function [v,a]=absorbance(v,L,v0,pcoeff,pshift,sigma,E,n,M,T,P,ppm,gammaLW)
% function [Pout,v]=icosspec(range,v0,pcoeff,sigma,M,T,P,conc)
% Input:
%     range: width of scan in wavenumbers 
%     Power: vector of power input to cavity
%     mirrorloss: mirror loss in ppm
%     L: cavity length
%     v0: vector of linecenters (cm-1)
%     pcoeff: Air-broadened halfwidth (HWHM) (cm-1/atm) @296K
%     pshift: Air-broadened pressure shift (cm-1/atm) @296K
%     sigma: vector of line strengths @296K (cm-1/(molecule cm-2)) 
%     E: Lower state energy (cm-1)
%     n: coeff. of temperature dependence
%     M: vector of atomic masses (amu)
%     T: temperature (K)
%     P: pressure (P)
%     ppm: mixing ratio (ppm)
%     gammaLW: Laser line width. To account for doppler broadening by laser line width. Can be 0 or left out
% Output:
%     Pout: simulated ICOS spectra 
%     v: vector of wavenumbers
%     
% ex. [Pout,v]=icosspec([1332,1334],[1332.3],[.0641],[5.76e-20],[16],[295],[20],[1.7]);
%     plot(v,Pout)

% v=[range(2):-1*(range(2)-range(1))/1999:range(1)];
%R = 1-mirrorloss*1e-6;
kb=1.380658e-23;        % Boltzmann's Constant in Joules/Kelvin
m=1.66e-27*M;           % convert from a.m.u. to kilograms
c=2.99792458e8;         % speed of light
c2 = 1.438789;          % second radiation constant c2 = hc/k (cm*K)
Tref=296;     % Reference temperature for Hitran database
%delta=10e-6;  %maximum fractional error in skew

%Temperature correction of line intensity
S=sigma.*(Tref/T).^(3/2).*exp(c2.*E*(T-Tref)/(Tref*T)).*(1-exp(-c2.*v0/T))./(1-exp(-c2.*v0/Tref));
%T and P correction of line halfwidth (gammaL)
gammaL=(Tref/T).^n.*pcoeff*P/760;           % Lorentzian HWHM
%Calculation of Doppler width
gammaD=v0.*sqrt(2*kb*T*log(2)./(m*c*c));   % Doppler HWHM
if nargin > 15
    gammaD=sqrt(gammaD.^2 + gammaLW^2);    % Doppler correction for laser line width
end
%Pressure-shift correction of line position
v0=v0 + pshift*P/760;
%Calculation of magnitude of voigt function
mag=ppm.*P./760*2.685e19./(T/273.15).*S*L./(sqrt(pi/log(2)).*gammaD);
voigtl=zeros(length(v0),length(v));
for i=1:length(v0)
    voigtl(i,:)=mag(i)*voigt3(v,v0(i),gammaD(i),gammaL(i))';
    if mod(i,10)==0
        disp(num2str(i));
    end
end
if length(v0)>1;
    a = sum(voigtl);
else
    a = voigtl;
end

return
