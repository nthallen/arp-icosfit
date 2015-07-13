function res2 = split_results(res1)
% First delete lines that have empty elements, and count lines
% that have non-scalar elements
flds = fields(res1);
nflds = length(flds);
nres1 = length(res1);
lengths = zeros(nflds, nres1);
for i=1:nflds
  fld = flds{i};
  lengths(i,:) = arrayfun(@(x) length(x{1}), {res1.(fld)});
end
hasres = min(lengths) > 0;
maxlen = max(lengths);
trows = sum(maxlen(hasres));
if trows == 0
  res2 = [];
else
  res2(trows) = res1(1);
  newrow = 1;
  for oldrow = 1:nres1
    if hasres(oldrow)
      if maxlen(oldrow) == 1
        res2(newrow) = res1(oldrow);
        newrow = newrow + 1;
      else
        mflds = lengths(:,oldrow)>1;
        if any(lengths(mflds,oldrow) ~= maxlen(oldrow))
          error('MATLAB:HUARP:MixedLengths', ...
            'Lengths do not match for oldrow %d', oldrow);
        end
        for mrow = 1:maxlen(oldrow)
          for i=1:nflds
            fld = flds{i};
            if mflds(i)
              res2(newrow).(fld) = res1(oldrow).(fld)(mrow);
            else
              res2(newrow).(fld) = res1(oldrow).(fld);
            end
          end
          newrow = newrow + 1;
        end
      end
    end
  end
end
