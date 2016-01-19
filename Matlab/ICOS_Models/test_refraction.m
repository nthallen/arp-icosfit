function test_refraction
% Sytematically look at refraction to test mathematical model
% I specifically want to look at skew and divergence, so I will
% assume Incident point I=(x,y,z) has z==0 OK to assume x==0 as well
% Also want axially symmetric optic, so the normal vector N=(X,Y,Z) has
% Z==0.
%
% To test then, we can select range of values for:
%   y, d, s, Y and calculate the rest.
%
% My goal is to figure out exactly how the approximations differ
% from the exact calculations and validate the equations I have worked
% through.
%
% Di = (x',y',z') is the unit incident direction.
%   si = incident skew = z'/x';
%   di = incident divergence = y'/x';
% N is unit normal vector, pointing forward
%   th_i is the incident angle from the normal
% ci = dot(Di,N) = cos(th_i) > 0
% Ni = ci*N = normal component of incident direction
% si = Di-Ni = component of incident direction perpendicular to normal
%      length of si = sin(th_i)
% ni = index of refraction on the incident side
% nt = index of refraction on the transmition side
% st = si*ni/nt = component of trasmitted direction perpendicular to normal
% ct = sqrt(1-len(st)^2) = length of normal component of transmitted
%       direction
% Nt = ct*N = normal component of trasmitted direction
% Dt = Nt+st = transmitted direction
  dsrange = 0.25;
  yrange = 0.2;
  ni = 1;
  nt = 2.4361;
  dv = linspace(0,dsrange,11);
  sv = linspace(0,dsrange,11);
  [dm,sm] = meshgrid(dv,sv);
  d = dm(:);
  s = sm(:);
  Di = Di_from_d_s(d,s);
  yv = linspace(0,yrange,11);
  for i = 1:length(yv)
    Y = yv(i) * ones(length(d),1);
    N = N_from_Y(Y);
    Dt = transmission(Di, N, ni, nt);
    % Verify that Dt are unit vectors
    lDt = std(sqrt(sum(Dt.^2,2))-1);
    if lDt > 1e-6
      warning('ARP:Normalize','Dt not unitized: lDt = %f at Y=%.3f', ...
        lDt, yv(i));
    end
    Dt2 = transmission2(Di, N, ni, nt);
    EDt2 = Dt2 - Dt;
    lDt2 = std(sqrt(sum(EDt2.^2,2)));
    if lDt2 > 1e-6
      warning('ARP:Normalize','Dt2 differs from Dt: lDt2 = %f at Y=%.3f', ...
        lDt2, yv(i));
    end
    % Assess skew approximation
    Sapprox = (ni/nt)*Di(:,3)./Di(:,1);
    ci = sum(Di.*N,2);
    Sexact = Di(:,3)./((sqrt((nt/ni)^2-1+ci.^2)-ci).*N(:,1) + Di(:,1));
    dS = (Sapprox-Sexact)./Sexact;
    
    % scatter(d,s,[],dS);
    
    dSm = reshape(dS,size(dm));
    mesh(dm,sm,dSm);
    xlabel('divergence');
    ylabel('skew');
    colorbar;
%    plot(d.*s,dS-.5*d.*s,'.');
    title(sprintf('Y = %.3f', yv(i)));
    shg;
    pause;
  end
end

function Di = Di_from_d_s(d, s)
  Di = [ones(size(d)),d,s];
  lDi = sqrt(sum(Di.^2,2));
  Di = Di./(lDi*ones(1,3));
end

function N = N_from_Y(Y)
  N = [sqrt(1-Y.^2),Y,zeros(size(Y))];
end

function Dt = transmission(Di, N, ni, nt)
  ci = sum(Di.*N,2);
  Ni = (ci*ones(1,3)).*N;
  Si = Di - Ni;
  St = Si*ni/nt;
  ct = sqrt(1-sum(St.^2,2));
  Dt = (ct*ones(1,3)).*N + St;
end
  
function Dt2 = transmission2(Di, N, ni, nt)
  ci = sum(Di.*N,2);
  nct = sqrt(nt^2/ni^2-1+ci.^2);
  Dt2 = (ni/nt)*[
    (nct-ci).*N(:,1) + Di(:,1), ...
    (nct-ci).*N(:,2) + Di(:,2), ...
    Di(:,3)
    ];
end
