function V = fitlin( v, n );
% V = fitlin( v , n );
%  v  exponential decay plus offset
%  n  number of samples to offset for the fit (defaults to 1)
% V = [ bminus b bplus s2 a ];
%         b = V(2);
%         a = V(5);
%         tau = n*dt/log(b);
%         z = a/(1-b);

if size(v,2) ~= 1
  error(['Input vector v must be a column vector']);
end
npts = length(v);

if npts < 1+n
  error(['Input vector not long enough']);
end
y = v(1:(npts-n));
x = v((1+n):npts);
V = linfit3(x,y);
