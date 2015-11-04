function dN = evaluate_interleave(m,k)
%% scratch_151021.m. Look at interleave to evaluate tolerances
% Back to the Future Day
% Interleave is characterized by:
%  m = number of spots
%  k = interleave factor
%  where gcd(m,k) == 1
dN = [];
k_inv = find(mod((1:m-1)*k,m)==1);
if isempty(k_inv)
  return;
end
% N = zeros(1,m+1);
% for i=0:m-1
%   Ni = mod(i*k,m);
%   N(Ni+1) = i;
% end
% N(m+1) = 0;
% dN = minmax(diff(N));
dN = minmax(diff(mod((0:m)*k_inv,m)));
