function check_params(sol_num, res)
% check_params(sol_num, res)
% Checks basic parameters against the ICOS equations
n = 2.4361;
run_check(sol_num, 'ICOS eq 1', check_eq1(res.h1, res.r2, res.R2, res.L));
run_check(sol_num, 'ICOS eq 2', check_eq1(res.h2, res.r1, res.R1, res.L));
run_check(sol_num, 'ICOS eq 3', check_eq3(res.r1, res.r2, res.h1, res.h2));
run_check(sol_num, 'ICOS RRL', check_rrl(res.R1, res.R2, res.L));
run_check(sol_num, 'ICOS M1 rhw', check_rhw(res.r1, res.h1, res.w1));
run_check(sol_num, 'ICOS M2 rhw', check_rhw(res.r2, res.h2, res.w2));
run_check(sol_num, 'ICOS s1', check_s(res.s1, res.w2, res.L, res.r2, res.R1, res.R2));
run_check(sol_num, 'ICOS s2', check_s(res.s2, res.w1, res.L, res.r1, res.R2, res.R1));
run_check(sol_num, 'ICOS d1', check_d(res.d1, res.r1, res.R1));
run_check(sol_num, 'ICOS d2', check_d(res.d2, res.r2, res.R2));
run_check(sol_num, 'ICOS Mirror 1', check_M1(res.r1, res.R1, res.Rr2, res.RR2, n));
run_check(sol_num, 'RIM eq 1', check_eq1(res.Rh1, res.Rr2, res.RR2, res.RL));
run_check(sol_num, 'RIM eq 2', check_eq1(res.Rh2, res.Rr1, res.RR1, res.RL));
run_check(sol_num, 'RIM eq 3', check_eq3(res.Rr1, res.Rr2, res.Rh1, res.Rh2));
run_check(sol_num, 'RIM RRL', check_rrl(res.RR1, res.RR2, res.RL));
run_check(sol_num, 'RIM M1 rhw', check_rhw(res.Rr1, res.Rh1, res.Rw1));
run_check(sol_num, 'RIM M2 rhw', check_rhw(res.Rr2, res.Rh2, res.Rw2));
run_check(sol_num, 'RIM s1', check_s(res.Rs1, res.Rw2, res.RL, res.Rr2, res.RR1, res.RR2));
run_check(sol_num, 'RIM s2', check_s(res.Rs2, res.Rw1, res.RL, res.Rr1, res.RR2, res.RR1));
run_check(sol_num, 'RIM d1', check_d(res.Rd1, res.Rr1, res.RR1));
run_check(sol_num, 'RIM d2', check_d(res.Rd2, res.Rr2, res.RR2));

function run_check(sol_num, test_desc, msgs)
for i=1:length(msgs)
  fprintf(1,'%d: %s: %s\n', sol_num, test_desc, msgs{i});
end

function msgs = check_eq1(h1, r2, R2, L)
% eq1 says h1/r2 = (R2-L)/R2. If R2-L==0
msgs = {};
Rrat = (R2-L)/R2;
if abs(Rrat) < .02
  hrat = h1/r2;
  if hrat == Rrat
    if Rrat == 0
      msgs{end+1} = 'R == L';
    else
      msgs{end+1} = 'R is very close to L';
    end
  elseif Rrat == 0
    msgs{end+1} = sprintf('h non-zero (%d) although R==L', h1);
  else
    msgs{end+1} = 'R is very close to L';
    msgs = compare_equal('h1',h1,'r2*(R2-L)/R2',Rrat*r2,1e-3, msgs);
  end
else
  msgs = compare_equal('r2',r2,'h1*R2/(R2-L)',h1/Rrat,1e-3, msgs);
end

function msgs = check_eq3(r1, r2, h1, h2)
msgs = {};
p = [r1, r2, h1, h2];
if any(p < 0)
  msgs{end+1} = 'Not all parameters are non-negative';
end
mp = min(abs(p));
Mp = max(abs(p));
if mp/Mp < .01
  msgs{end+1} = sprintf('Nearly degenerate solution: %.2%%', 100*mp/Mp);
  % should check 
else
  msgs = compare_equal('r1/r2', r1/r2, 'h1/h2', h1/h2, 1e-4, msgs);
end

function msgs = check_rhw(r, h, w)
msgs = {};
if h > r
  msgs{end+1} = 'h > r';
end
if w > r
  msgs{end+1} = 'w > r';
end
msgs = compare_equal('r^2', r^2, 'h^2+w^2', h^2+w^2, 1e-3, msgs);

function msgs = check_s(s1, w2, L, r2, R1, R2)
msgs = compare_equal('s1', s1, 'w2/L', w2/L, 1e-3, {});
if (R1+R2-L)/(R1*R2) <= 0
  msgs{end+1} = '(R1+R2-L)/(R1*R2) < 0';
else
  w2a2 = r2^2*L*(R1+R2-L)/(R1*R2);
  msgs = compare_equal('w2^2', w2^2, 'f(r2,R1,R2,L)', w2a2, 1e-3, msgs);
end

function msgs = check_d(d, r, R)
msgs = compare_equal('d', d, 'r/R', r/R, 1e-3, {});

function msgs = check_M1(r1, R1, Rr2, RR2, n)
msgs = compare_equal('r1', r1, 'Rr2', Rr2, 0, {});
msgs = compare_equal('RR2', RR2, '-R1/n', -R1/n, 0, msgs);

function msgs = check_rrl(R1, R2, L)
msgs = {};
if (R1+R2-L)/(R1*R2) <= 0
  msgs{end+1} = '(R1+R2-L)/(R1*R2) <= 0';
end

function msgs = compare_equal(l1, v1, l2, v2, threshold, msgs)
rdiff = (v1-v2)/v1;
if abs(rdiff) > threshold
  msgs{end+1} = sprintf('%s differs from %s by %.2f%%', ...
    l2, l1, 100*rdiff);
end
