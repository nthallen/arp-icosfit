function T = cpci14time( cpci14 );
% T = cpci14time( cpci14 );
PT = load_crmat('PT');
v = find(diff(PT.CPCI14)>0)+1;
T = interp1(PT.CPCI14(v), PT.TPT(v), cpci14, 'linear', 'extrap');
