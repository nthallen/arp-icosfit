function res = collect_results(varargin)
% res = collect_results(params);
% ICOS_search objects named 'IS'.
% params is a struct. The fields are the fields from the IS.res2
% struct. The values, if present, identify which range of values
% to include in the results collection.
%%
files = dir('IS_*.mat'); % { 'IS_w30_L50.mat', 'IS_w35_L50.mat' };
files = { files.name };
isIB = regexp(files,'(\dx\d)|(_save)');
cellisempty = @(x) isempty(x{1});
isIB2 = arrayfun(cellisempty, isIB);
files = files(isIB2)';
%%
params.RR1 = [];
params.RD1 = [];
params.Rw1 = [];
params.RL = [];
params.r1 = [];
params.L = [];
params.R2 = [];
params.D2 = [];
params.overlap = [];
params.theta = [];
params.r_d = [];
params.NH = [];
params.max_pwr = [];
%%
for i=1:2:length(varargin)-1
  fld= varargin{i};
  params.(fld) = varargin{i+1};
end

pflds = fieldnames(params);
nfiles = length(files);
col(nfiles).ri = [];
for i=1:length(files)
  IS = load(files{i});
  IS = IS.IS;
  col(i).IS = IS;
  col(i).ri = zeros(1,length(IS.res2)); % row to facilitate aggregation
  for j=1:length(IS.res2)
    OK = 1;
    for k=1:length(pflds)
      if OK
        fld = pflds{k};
        if ~isfield(IS.res2(j), fld) || isempty(IS.res2(j).(fld)) % must be non-empty, so analyzed
          OK = 0;
        elseif ~isempty(params.(fld)) % must be non-empty, so analyzed
          cons = params.(fld);
          if length(cons) == 1
            if IS.res2(j).(fld) ~= cons
              OK = 0;
            end
          elseif length(cons) == 2
            if (IS.res2(j).(fld) < params.(fld)(1) || ...
                IS.res2(j).(fld) > params.(fld)(2))
              OK = 0;
            end
          else
            error('MATLAB:HUARP:BadConstraint', ...
              'Constraint for parameter %s has unsupported length %d', ...
              fld, length(cons));
          end
        end
      end
    end
    col(i).ri(j) = OK;
  end
end
% Now add up the number of results, allocate the array and fill it.
%%
nres = sum([col.ri]);
sargs = cell(1,2*(length(pflds)+2));
sargs{1} = 'mnc';
sargs{3} = 'index';
sargs(5:2:2*(length(pflds)+2)) = pflds;
res(nres) = struct(sargs{:});
%%
result = 1;
for i=1:length(files)
  ri = find(col(i).ri);
  if iscolumn(ri)
    ri = ri';
  end
  for j=ri
    res(result).mnc = col(i).IS.ISopt.mnc;
    res(result).index = j;
    for k=1:length(pflds)
      fld = pflds{k};
      res(result).(fld) = col(i).IS.res2(j).(fld);
    end
    result = result+1;
  end
end
%%
for i=1:length(res)
  %%
  pat1 = sprintf('IS_%s.%d_*x*.mat', res(i).mnc, res(i).index);
  IBs = dir(pat1);
  IBs = { IBs.name };
  IBsana = regexp(IBs,'\d+x\d+\.mat');
  IBs = IBs(~isempty(IBsana));
  if isempty(IBs)
    fprintf(1,'No ICOS_beam analysis found for IS_%s.%d\n', ...
      res(i).mnc, res(i).index);
  else
    if length(IBs) > 1
      fprintf(1,'More than one ICOS_beam analysis found for IS_%s.%d\n', ...
        res(i).mnc, res(i).index);
    end
    load(IBs{1});
    Pwr = IB.PowerSummary;
    res(i).eNH = res(i).NH * (1-Pwr.H_loss/100);
    res(i).H_loss = Pwr.H_loss;
    res(i).I_loss = Pwr.I_loss;
    res(i).F_loss = Pwr.F_loss;
    res(i).D_loss = Pwr.D_loss;
  end
end
