function [ emdl, J] = etln_evalJ(Y,x)
% [ emdl, J ] = etln_evalJ(Y,x);
% Models the etalon airy function for a given sample.
% Y is the parameter vector
%   The first 5 or 7 parameters are those for etln_eval() and
%   determine the fringenumber at each sample
%   The next 4 parameters are the power curve 3rd order fit
%   The last parameter is the etalon Finesse
% rx is the sample number vector
x2 = x.*x;
if length(Y) == 5
  cse1 = exp(-x/Y(5));
  cse2 = Y(4)*cse1;
  emdl = Y(1)+Y(2)*x+Y(3)*x2+cse2;
  if nargout > 1
    J = zeros(length(x),length(Y));
    J(:,1) = 1;
    J(:,2) = x;
    J(:,3) = x2;
    J(:,4) = cse1;
    J(:,5) = cse2.*x/(Y(5)*Y(5));
  end
elseif length(Y) == 7
  cse1 = exp(-x/Y(5));
  cse2 = Y(4)*cse1;
  cse3 = exp(-x/Y(7));
  cse4 = Y(6)*cse3;
  emdl = Y(1)+Y(2)*x+Y(3)*x2+cse2+cse4;
  if nargout > 1
    J = zeros(length(x),length(Y));
    J(:,1) = 1;
    J(:,2) = x;
    J(:,3) = x2;
    J(:,4) = cse1;
    J(:,5) = cse2.*x/(Y(5)*Y(5));
    J(:,6) = cse3;
    J(:,7) = cse4.*x/(Y(7)*Y(7));
  end
elseif length(Y) == 10 % 5 for fringe number
  x3 = x2.*x;
  P = Y(6)*x3 + Y(7)*x2 + Y(8)*x + Y(9);
  if nargout > 1
    [ f, Jf ] = etln_evalJ( Y(1:5),x);
  else
    f = etln_evalJ( Y(1:5),x );
  end
  cse5 = pi*f;
  cse6 = sin(cse5);
  cse7 = cse6.*cse6;
  D = 1 + Y(10)*cse7;
  emdl = P./D;
  % plot(x,emdl); shg; pause;
  if nargout > 1
    J = zeros(length(x),length(Y));
    J(:,1:5) = - emdl*2*pi*Y(10).*cse6.*cos(cse5)*ones(1,5) .* Jf;
    J(:,6)  = x3;
    J(:,7)  = x2;
    J(:,8) = x;
    J(:,9) = 1;
    J(:,10) = - emdl.*cse7;
    J = J./(D*ones(1,10));
  end
else % assume length 12
  x3 = x2.*x;
  P = Y(8)*x3 + Y(9)*x2 + Y(10)*x + Y(11);
  if nargout > 1
    [ f, Jf ] = etln_evalJ( Y(1:7),x);
  else
    f = etln_evalJ( Y(1:7),x );
  end
  cse5 = pi*f;
  cse6 = sin(cse5);
  cse7 = cse6.*cse6;
  D = 1 + Y(12)*cse7;
  emdl = P./D;
  % plot(x,emdl); shg; pause;
  if nargout > 1
    J = zeros(length(x),length(Y));
    J(:,1:7) = - emdl*2*pi*Y(12).*cse6.*cos(cse5)*ones(1,7) .* Jf;
    J(:,8)  = x3;
    J(:,9)  = x2;
    J(:,10) = x;
    J(:,11) = 1;
    J(:,12) = - emdl.*cse7;
    J = J./(D*ones(1,12));
  end
end
