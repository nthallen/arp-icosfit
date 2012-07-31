function pick_regions
% pick_regions
% Strategy:
%  1: identify contiguous regions in each waveform
%  2: within ringdown regions, restrict to regions where
%     pinch valve is open
%  3: within icos regions, identify two baseline regions and a sample
%     region
line_obj = fitline('load');
waves = load_waves;
D = load('creng_1.mat');
x = [1:length(D.CPCI14)]';
v = find(diff(D.CPCI14)>0)+1; % Index of CPCI numbers
wn = D.QCLI_Wave(v);
vv = [ 0; find(diff(wn)) ] + 1; % Index in v of new waveforms
wnum = wn(vv);
wnstart = D.CPCI14(v(vv)); % CPCI of start
wnend = [ D.CPCI14(v(vv(2:end)))-1; max(D.CPCI14) ]; % CPCI of end
xstart = v(vv);
xend = v([vv(2:end); end])-1;

% find ringdown regions
% rdreg = find(~[ waves(wnum).ISICOS ]);
% for i=1:length(rdreg)
%   line_obj = addregion( line_obj, 'RD', i, [wnstart(rdreg(i)) wnend(rdreg(i))]);
% end

% identify Sample Regions
icosreg = find([waves(wnum).ISICOS]);
samplenum = 0;
for i=1:length(icosreg)
  % find index into the D arrays where we have new cpci numbers
  % within this region
  xreg = x(v(x(v) >= xstart(icosreg(i)) & x(v) <= xend(icosreg(i))));
  % Now find the places where the pinch valve is open
  % And CellP > 20
  % ibit, samplestart and sampleend are all referenced to
  % xreg
  ibit = bitand(D.PVStat(xreg),1) == 0;
  [samplestart,sampleend] = select_reg( ibit, 10 );
  for j=1:length(samplestart)
    samplereg = xreg(samplestart(j):sampleend(j));
    vdP = abs(diff(D.CellP(samplereg))) < .4;
    [Pstart,Pend] = select_reg(vdP, 10);
    if ~isempty(Pstart)
      samplenum = samplenum + 1;
      samplePreg = samplereg(Pstart(1):end);
      cpcistart = D.CPCI14(samplePreg(1));
      cpciend = D.CPCI14(samplePreg(end));
      line_obj = addregion( line_obj, 'R', samplenum, [ cpcistart cpciend ] );
      vdT = [ 1; diff(D.Tcreng_1(samplePreg)) < 10 ];
      [qclistart,qcliend] = select_reg(vdT,10);
      if length(qclistart) > 1 || qclistart(1) > 1
        for j1 = 1:length(qclistart)
          Tqcli = D.Tcreng_1(samplePreg(qclistart(j1):qcliend(j1)));
          if qclistart(j1) > 1
            iQ = min(find(Tqcli > Tqcli(1)+60)); % QCLI Warmup
            if isempty(iQ)
              qclistart(j1) = qcliend(j1);
            else
              qclistart(j1) = qclistart(j1) + iQ - 1;
            end
          end
          if qclistart(j1) < qcliend(j1)
            line_obj = addregion( line_obj, 'R', samplenum, ...
              D.CPCI14(samplePreg([qclistart(j1) qcliend(j1)])), j1 );
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
%     cpciBstart = D.CPCI14(xBstart);
%     cpciBend = D.CPCI14(xBend);
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
%     cpciRstart = D.CPCI14(xRstart);
%     cpciRend = D.CPCI14(xRend);
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
  lo_in.Regions(end+1) = struct('name', name, 'cpci', range );
  [ names, ndx ] = sort( {lo_in.Regions.name} );
  lo_in.Regions = lo_in.Regions(ndx);
  lo_in.CurRegion = find(strcmp({lo_in.Regions.name}, 'all'));
  if length(lo_in.CurRegion) ~= 1
    lo_in.CurRegion = 1;
  end
end
lo_out = lo_in;
return
