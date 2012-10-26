function T = scantime( scannum );
% T = scantime( cpci14 );
PT = load_mat_files('PT');
v = find(diff(PT.ScanNum)>0)+1;
T = interp1(PT.ScanNum(v), PT.TPT(v), scannum, 'linear', 'extrap');
