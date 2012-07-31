function writebin(ofile, fe);
% writebin( ofile, fe );
% Writes the matrix fe to the file ofile in the
% standard binary format for Cavity Ringdown.
% i 4 nrows
% i 4 ncols
% REPEAT ncols {
%   REPEAT nrows {
%     r 4 value
%   }
% }
fid = fopen(ofile, 'w');
if fid > 0
  fwrite(fid, size(fe), 'integer*4' );
  fwrite(fid, fe, 'real*4' );
  fclose(fid);
else
  error(['Error writing file: ' ofile ]);
end
