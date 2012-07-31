function fe = loadbin( ifile );
% fe = loadbin( ifile );
% returns an empty matrix if the file cannot be read.
fid = fopen(ifile, 'r');
if fid > 0
  try
    [ dim, count ] = fread(fid,2,'uint32','l');
    fe = fread(fid,dim','float32','l');
  catch
    fe = [];
  end
  fclose(fid);
else
  fe = [];
end
