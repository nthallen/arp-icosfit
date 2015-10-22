%%
% Demonstrate different interleave patterns
C = 3000;
W = .4;
R = 5;
rd = .1;
th = 15;
L = 90;
Nf = C/(2*L);
I = 1;
N = ceil(Nf); % This spot can overlap
wp = W/2;
for NE = 4:-1:0
  phi1 = pi/(N+NE); % half angle of single spot spacing
  phiI = I*phi1;
  w = wp*sin(phiI)/sin(phi1);
  r = wp/sin(phi1);
  draw_spots(R, r, phiI*2, N);
  pause;
end
I = 2;
NI = 2*ceil((Nf-1)/2)+1; % smallest number of spots to fit N spots
for NE = (N-1):-2:0
  phi1 = pi/(N+NE); % half angle of single spot spacing
  phiI = I*phi1;
  w = wp*sin(phiI)/sin(phi1);
  r = wp/sin(phi1);
  draw_spots(R, r, phiI*2, N);
  pause;
end
%%
% Another approach
I = 2;
NI = I*ceil((Nf-1)/I)+1; % smallest number of spots to fit N spots
for NO = 1:NI
  if gcd(NO,I) == 1
    phi1 = pi/(I*N-NO); % half angle of single spot spacing
    phiI = I*phi1;
    w = wp*sin(phiI)/sin(phi1);
    r = wp/sin(phi1);
    draw_spots(R, r, phiI*2, N);
    title(sprintf('I:%d N:%d NO:%d', I, N, NO));
    pause;
  end
end
%%
I = 4;
NI = I*ceil((Nf-1)/I)+1; % smallest number of spots to fit N spots
for NO = 1:N*(I-1)
  if gcd(NO,I) == 1
    phi1 = pi/(I*N-NO); % half angle of single spot spacing
    phiI = I*phi1;
    w = wp*sin(phiI)/sin(phi1);
    r = wp/sin(phi1);
    draw_spots(R, r, phiI*2, N);
    title(sprintf('I:%d N:%d NO:%d', I, N, NO));
    pause;
  end
end
