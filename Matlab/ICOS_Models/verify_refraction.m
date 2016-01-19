function verify_refraction
  R = 75; % cm
  n1 = 1;
  n2 = 2.4361; % ZiSn
  r_max = 4; % cm
  theta_max = 2*pi/180; % radians
  d_max = tan(theta_max);
  N = 100;
  I = rand(N,4)*diag([r_max,r_max,d_max,d_max]); % Incident beam
  M = [
    1 0 (1-n1/n2)/R 0
    0 1 0 (1-n1/n2)/R
    0 0 n1/n2 0
    0 0 0 n1/n2];
  R1 = I*M; % 4-component vector
  % Now do it the classical way
  DV = [ones(N,1),I(:,[3,4])]; % direction vector of incident
  IV = [R-sqrt(R^2-I(:,1).^2 - I(:,2).^2),I(:,[1 2])]; % Intersection point
  NV = (IV - ones(N,1)*[R,0,0])/R;
  NRV = sum(DV.*NV,2)*ones(1,3).*NV; % Normal component of DV
  PV = DV - NRV; % perpendicular component
  R2 = NRV + PV*n1/n2;
  
  R1V = [ones(N,1),R1(:,[3,4])];
  EV = R1V-R2;
  MEV = sqrt(sum(EV.^2,2));
  r = sqrt(sum(I(:,[1,2]).^2,2));
  dr = sqrt(sum(I(:,[3,4]).^2,2));
  figure;
  plot(r,MEV,'*');
  xlabel('r, cm');
  ylabel('Error, cm');
  figure;
  plot(r,atand(MEV),'*');
  ylabel('Error degrees');
  xlabel('r, cm');
  figure;
  plot(dr,atand(MEV),'*');
  xlabel('dr, cm');
  ylabel('Error degrees');
  
  d1 = (R1(:,1).*R1(:,3)+R1(:,2).*R1(:,4))./sqrt(sum(I(:,[1,2]).^2,2));
  s1 = (R1(:,1).*R1(:,4)-R1(:,2).*R1(:,3))./sqrt(sum(I(:,[1,2]).^2,2));
  
  d2 = (I(:,1).*R2(:,2)+I(:,2).*R2(:,3))./sqrt(sum(I(:,[1,2]).^2,2));
  s2 = (I(:,1).*R2(:,3)-I(:,2).*R2(:,2))./sqrt(sum(I(:,[1,2]).^2,2));
  
  figure;
  plot(r,d1-d2,'*');
  xlabel('r,cm');
  ylabel('Divergence error');
  
  figure;
  plot(r,s1-s2,'*');
  xlabel('r,cm');
  ylabel('skew error');
  
  
