function pick_regions
% pick_regions
% Strategy:
%  1: identify contiguous regions in each waveform
%  2: within ringdown regions, restrict to regions where
%     pinch valve is open
%  3: within icos regions, identify two baseline regions and a sample
%     region
% 
Algo='A';
cfg = load_ICOSfit_cfg;
line_obj = fitline('load');
waves = load_waves;
D = ne_load('HCIeng_1', 'HCI_Data_Dir');
D.CPrs_Reg_DS = bitand(D.DS84C,32);
D.CGas_Vlv_DS = bitand(D.DS84C,64);
D.MGas_Vlv_DS = bitand(D.DS868,128);
D.MPrs_Reg_DS = bitand(D.DS868,64);
SSP_Num = eval(['D.' cfg.ScanDir '_Num']);
x = [1:length(SSP_Num)]';
v = find(diff(SSP_Num)>0)+1; % Index of SSP numbers
QCLI_Wave = eval(['D.' cfg.WavesFile(1:6) '_Wave']);
wn = QCLI_Wave(v);
vv = [ 0; find(diff(wn)) ] + 1; % Index in v of new waveforms
wnum = wn(vv);
wnstart = SSP_Num(v(vv)); % SSP of start
wnend = [ SSP_Num(v(vv(2:end)))-1; max(SSP_Num) ]; % SSP of end
xstart = v(vv);
xend = v([vv(2:end); end])-1;

% find ringdown regions
% rdreg = find(~[ waves(wnum).ISICOS ]);
% for i=1:length(rdreg)
%   line_obj = addregion( line_obj, 'RD', i, [wnstart(rdreg(i)) wnend(rdreg(i))]);
% end

% identify Sample Regions
icosreg = find([waves(wnum+1).ISICOS]);
samplenum = 0;
Axis=cfg.ScanDir(5);
if strcmp(Axis,'I')
    Axis='M';
end
if strcmp(Algo,'A')
    Axis='C';
end
for i=1:length(icosreg)
  % find index into the D arrays where we have new SSP numbers
  % within this region
  xreg = x(v(x(v) >= xstart(icosreg(i)) & x(v) <= xend(icosreg(i))));
  % Now find the places where the solenoid valve is open
  % And CellP > 30
  % ibit, samplestart and sampleend are all referenced to
  % xreg
  
  PRStat = eval(['D.' Axis 'Prs_Reg_DS']);
  CellP=eval(['D.' Axis 'CelP']);
  ibit = bitand(PRStat(xreg),max(PRStat)) ~= 0;
  [samplestart,sampleend] = select_reg( ibit, 150 );
  for j=1:length(samplestart)
    samplereg = xreg(samplestart(j):sampleend(j));
    vdP = abs(diff(CellP(samplereg))) < .4;
    [Pstart,Pend] = select_reg(vdP, 10);
    if ~isempty(Pstart)
      samplenum = samplenum + 1;
      samplePreg = samplereg(Pstart(1):end);
      sspstart = SSP_Num(samplePreg(1));
      sspend = SSP_Num(samplePreg(end));
      line_obj = addregion( line_obj, 'R', samplenum, [ sspstart sspend ] );
      vdT = [ 1; diff(D.THCIeng_1(samplePreg)) < 10 ];
      [qclistart,qcliend] = select_reg(vdT,10);
      if length(qclistart) > 1 || qclistart(1) > 1
        for j1 = 1:length(qclistart)
          Tqcli = D.THCIeng_1(samplePreg(qclistart(j1):qcliend(j1)));
          if qclistart(j1) > 1
            iQ = find(Tqcli > Tqcli(1)+60, 1 ); % QCLI Warmup
            if isempty(iQ)
              qclistart(j1) = qcliend(j1);
            else
              qclistart(j1) = qclistart(j1) + iQ - 1;
            end
          end
          if qclistart(j1) < qcliend(j1)
            line_obj = addregion( line_obj, 'R', samplenum, ...
              SSP_Num(samplePreg([qclistart(j1) qcliend(j1)])), j1 );
          end
        end
      end
    end
  end
  PRStat = eval(['D.' Axis 'Gas_Vlv_DS']);
  ibit = bitand(PRStat(xreg),max(PRStat)) ~= 0;
  [samplestart,sampleend] = select_reg( ibit, 30 );
  for j=1:length(samplestart)
    samplereg = xreg(samplestart(j):sampleend(j));
    vdP = abs(diff(CellP(samplereg))) < 3;
    [Pstart,Pend] = select_reg(vdP, 10);
    if ~isempty(Pstart)
      %samplenum = samplenum + 1;
      samplePreg = samplereg(Pstart(1):end);
      sspstart = SSP_Num(samplePreg(1));
      sspend = SSP_Num(samplePreg(end));
      line_obj = addregion( line_obj, 'C', samplenum, [ sspstart sspend ] );
      vdT = [ 1; diff(D.THCIeng_1(samplePreg)) < 10 ];
      [qclistart,qcliend] = select_reg(vdT,10);
      if length(qclistart) > 1 || qclistart(1) > 1
        for j1 = 1:length(qclistart)
          Tqcli = D.THCIeng_1(samplePreg(qclistart(j1):qcliend(j1)));
          if qclistart(j1) > 1
            iQ = find(Tqcli > Tqcli(1)+60, 1 ); % QCLI Warmup
            if isempty(iQ)
              qclistart(j1) = qcliend(j1);
            else
              qclistart(j1) = qclistart(j1) + iQ - 1;
            end
          end
          if qclistart(j1) < qcliend(j1)
            line_obj = addregion( line_obj, 'C', samplenum, ...
              SSP_Num(samplePreg([qclistart(j1) qcliend(j1)])), j1 );
          end
        end
      end
    end
  end
end
save fitline.mat line_obj
return

% RegNum = 0;
% xcalstart = find(diff(D.CalSt>0)>0)+1;
% xcalend = find(diff(D.CalSt>0)<0);
% PVStep = D.PVStep;
% vPV = PVStep > 60000;
% PVStep(vPV) = PVStep(vPV) - 65536;
% xPstart = find(diff(PVStep>0)>0)+25;
% xPend = find(diff(PVStep>0)<0)-2;
% Regions = struct('name', 'all', 'cpci', [] );
% for i = 1:length(wnum)
%   if (waves(wnum(i)).ISICOS)
%     % pick two Bs and one R
%     xBstart = xcalstart(find(xcalstart > xstart(i) & xcalstart < xend(i)))+5;
%     xBend = xcalend(find(xcalend > xstart(i) & xcalend < xend(i)));
%     if length(xBstart) == 0
%       xBstart = xstart(i);
%     end
%     if length(xBend) == 0
%       xBend = xend(i);
%     end
%     if xBend(1) < xBstart(1)
%       xBstart = [ xstart(i); xBstart ];
%     end
%     if length(xBstart) > length(xBend)
%       xBend = [ xBend; xend(i) ];
%     end
%     if length(xBstart) ~= length(xBend)
%       error('Lengths do not match');
%     end
%     if any(xBstart > xBend)
%       error('Bad ranges');
%     end
%     cpciBstart = SSP_Num(xBstart);
%     cpciBend = SSP_Num(xBend);
%     for j=1:length(cpciBstart)
%       regname = sprintf('B%d%c', RegNum, 'a'+j-1);
%       range = [cpciBstart(j) cpciBend(j)];
%       Regions = [ Regions; struct( 'name', regname, 'cpci', range ) ];
%     end
%     xRstart = xPstart(find(xPstart > xstart(i) & xPstart < xend(i)));
%     xRend = xPend(find(xPend > xstart(i) & xPend < xend(i)));
%     if length(xRstart) == 0
%       xRstart = xstart(i);
%     end
%     if length(xRend) == 0
%       xRend = xend(i);
%     end
%     if xRend(1) < xRstart(1)
%       xRstart = [ xstart(i); xRstart ];
%     end
%     if length(xRstart) > length(xRend)
%       xRend = [ xRend; xend(i) ];
%     end
%     if length(xRstart) ~= length(xRend)
%       error('Lengths do not match');
%     end
%     if any(xRstart > xRend)
%       error('Bad ranges');
%     end
%     cpciRstart = SSP_Num(xRstart);
%     cpciRend = SSP_Num(xRend);
%     if length(cpciRstart) > 1
%       error('Not expecting two R regions');
%     end
%     for j=1:length(cpciRstart)
%       regname = sprintf('R%d', RegNum );
%       range = [cpciRstart(j) cpciRend(j)];
%       Regions = [ Regions; struct( 'name', regname, 'cpci', range ) ];
%     end
%   else
%     % pick an RD
%     RegNum = RegNum + 1;
%     
%   end
% end
% [Y,I] = sort({Regions.name});
% Regions = Regions(I);

function [ r_start, r_end ] = select_reg( v, minlen )
samplestart = find( diff([ 0; v ] ) > 0 );
sampleend = find( diff( [ v; 0 ] ) < 0 );
if length(samplestart) ~= length(sampleend)
  error('select_reg had some trouble');
end
len = sampleend-samplestart+1;
if any(len<=0)
  error('Came up with negative length in select_reg');
end
r_start = samplestart(len >= minlen);
r_end = sampleend(len >= minlen);
return

function lo_out = addregion( lo_in, base, num, range, suffix );
name = sprintf( '%s%d', base, num );
if nargin > 4
  name = [ name char('a'+suffix-1) ];
end
if diff(range) <= 0
  error('Region %s has negative duration', name);
end
if any(strcmp({lo_in.Regions.name}, name))
  fprintf( 1, 'Region %s [ %d %d ]: Not replaced\n', name, range );
else
  fprintf( 1, 'Region %s [ %d %d ]\n', name, range );
  lo_in.Regions(end+1) = struct('name', name, 'scan', range );
  [ names, ndx ] = sort( {lo_in.Regions.name} );
  lo_in.Regions = lo_in.Regions(ndx);
  lo_in.CurRegion = find(strcmp({lo_in.Regions.name}, 'all'));
  if length(lo_in.CurRegion) ~= 1
    lo_in.CurRegion = 1;
  end
end
lo_out = lo_in;
return
