function T = scantime( cpci14 );
% T = scantime( cpci14 );
PT = load_mat_files('PT');
v = find(diff(PT.CPCI14)>0)+1;
T = interp1(PT.CPCI14(v), PT.TPT(v), cpci14, 'linear', 'extrap');
