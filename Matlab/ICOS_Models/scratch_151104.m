%%
% Investigate different interleave patterns and evaluate for relevant
% design criteria:
%   r_min = r1_min (min pattern radius to avoid overlap)
%   r_max = r1_max (max pattern radius to achieve focus)
C = 3000;
B = .4;
rd = .1;
th = 15;
L = 35;
r_limit = (1.55*2.54 - B); % 2.5" diameter with beam margin
m_min = ceil(C/(2*L)); % This spot can overlap
m_max = floor(pi/asin(B/(2*r_limit)));
draw_all = false;
SolInc = 40;
RawSols = SolInc;
NSol = 0;
Summary(RawSols) = ...
  struct('RLmin',[],'r_max',[],'R1',[],'R2',[],'m',[],'k',[],'Phi',[], ...
    'Phi_n', [], 'Phi_p', [], 'Phi_N', [], 'Phi_P', [], 'Phi_OK', []);
%%
rdtanth = rd*tand(th);
figure;
for m=m_min:m_max
  for k = 1:floor(m/2);
    if gcd(m,k) ~= 1
      continue;
    end
    phi = pi/m;
    Phi = k*phi;
    R1pole = L/sin(Phi)^2;
    if R1pole > 2000
      R1 = interp1([0 2000],[L 2000],1:2000)';
      R2 = L*(R1-L)./(sin(Phi)^2*R1-L);
      RR = []; % L/(1+cos(Phi));
      Rmax = 2000;
    else
      R1a = interp1([0 1001],[L R1pole],1:1000)';
      R1b = interp1([0 1000],[R1pole 2000],1:1000)';
      R1 = [R1a;R1pole;R1b];
      R2 = [L*(R1a-L)./(sin(Phi)^2*R1a-L); NaN; L*(R1b-L)./(sin(Phi)^2*R1b-L)];
      RR = L./(1+cos(Phi)*(-1));
      Rmax = min(RR*1.1,2000);
      RR = RR(RR<2000);
    end
    R2(R2 > 2000 | R2 < -2000) = NaN;
    % fprintf(1,'Break at R1 = %.1f\n', L/sin(Phi)^2);
    r_min = B/(2*sin(pi/m));
    r_max = sqrt(rdtanth*sqrt(((R2-L).*R1.^2*L)./((R1-L).*(R1+R2-L))));
    r_max(~isnan(r_max)) = min(r_max(~isnan(r_max)),r_limit);
    rRR = interp1(R1,r_max,RR,'linear','extrap');
    % This (commented) calculation is wrong, but I haven't figured out why yet
    %   k2 = ((R1-L).*(R1.*R2-L))./((R2-L).*R1.^2.*L);
    %   k2(k2<0) = NaN;
    %   s1 = r_max.*sqrt(k2);
    r2 = abs(r_max.*cos(Phi).*R2./(R2-L));
    w2 = abs(r_max*sin(Phi)*cos(Phi).*R2./(R2-L));
    s1 = w2/L;
    
    Rw1 = 1.1*B/2;
    RL = abs(Rw1./s1);
    rRL = interp1(R1,RL,RR,'linear','extrap');
    vok = r_max > r_min & r_max >= r2;
    vnok = ~vok;
    RLok = RL;
    RLok(vnok) = NaN;
    RLnok = RL;
    RLnok(vok) = NaN;
    
    NSol = NSol+1;
    if NSol > RawSols
      RawSols = RawSols + SolInc;
      Summary(RawSols).RLmin = [];
    end
    
    dN = evaluate_interleave(m,k);
    
    Summary(NSol).RLmin = nanmin(RLok);
    mi = find(RLok == Summary(NSol).RLmin);
    phi_min = asin(B./r_max(mi));
    Phi_n = Phi + (k*phi_min-Phi)/dN(2);
    Phi_p = Phi + (k*phi_min-Phi)/dN(1);
    Summary(NSol).R1 = R1(mi);
    Summary(NSol).R2 = R2(mi);
    Summary(NSol).r_max = r_max(mi);
    Summary(NSol).m = m;
    Summary(NSol).k = k;
    Summary(NSol).Phi = Phi;
    Summary(NSol).Phi_n = Phi_n;
    Summary(NSol).Phi_p = Phi_p;
    SP.R1 = R1(mi);
    SP.R2 = R2(mi);
    SP.L = L;
    SP.Rw1 = 1.1*B/2;
    SP.RL = SP.Rw1/s1(mi);
    Res = exparam(SP);
    check_params(NSol, Res);
    
    if draw_all
      clf;
      ax = [nsubplot(3,1,1) nsubplot(3,1,2) nsubplot(3,1,3)];
      plot(ax(1),R1,r_max,'b',[R1(1) R1(end)],r_min*[1 1],'g',RR,rRR,'*r',...
        R1,r2);
      % xlabel(ax(1),'R1 cm');
      ylabel(ax(1),'r cm');
      %xlim(ax(1),[L/2 Rmax]);
      ylim(ax(1),[0 5]);
      set(ax(1),'XTickLabel',[],'YAxisLocation','Right');
      title(ax(1),sprintf('L=%d %d/%d', L, k, m));
      
      plot(ax(3),R1,R2);
      ylabel(ax(3),'R_2');
      set(ax(3),'YAxisLocation','Right');
      plot(ax(2),R1,RLok,'b',R1,RLnok,'r',RR,rRL,'*r');
      linkaxes(ax,'x');
      set(ax,'xlim',[L/2 Rmax]);
      ylim(ax(2),[0 min(L,nanmax(RL))]);
      xlabel(ax(end),'R1 cm');
      ylabel(ax(2),'RL cm');
      set(ax(2),'XTickLabel',[],'YAxisLocation','Left');
      % title(sprintf('RL: L=%d %d/%d', L, k, m));
      drawnow;
      pause;
    end
  end
end
%%
Summary = Summary(1:NSol);
%%
for i=1:length(Summary)
  R1 = Summary(i).R1 * (1+.01*[-1,1,-1,1]);
  R2 = Summary(i).R2 * (1+.01*[-1,-1,1,1]);
  Summary(i).Phi_OK = false;
  Phi_a = L.*(R1+R2-L)./(R1.*R2);
  if all(Phi_a >= 0 & Phi_a <= 1)
    Phi = asin(sqrt(Phi_a));
    Summary(i).Phi_N = Phi(1);
    Summary(i).Phi_P = Phi(2);
    Summary(i).Phi_OK = Phi(1)>= Summary(i).Phi_n & Phi(2) <= Summary(i).Phi_p;
  end
end
%%
clf;
ax = [ nsubplot(3,1,1), nsubplot(3,1,2), nsubplot(3,1,3)];
Phi = [Summary.Phi];
RLmin = [Summary.RLmin];
Phi_OK = [Summary.Phi_OK];
r_max = [Summary.r_max];
plot(ax(1),Phi(~Phi_OK), RLmin(~Phi_OK), 'r.',Phi(Phi_OK), RLmin(Phi_OK), 'b*');
plot(ax(2),Phi(~Phi_OK), r_max(~Phi_OK), 'r.',Phi(Phi_OK), r_max(Phi_OK), 'b*');
plot(ax(3),[Summary.Phi], [Summary.R1], '*', [Summary.Phi], [Summary.R2], '+');
title(ax(1),sprintf('L=%d', L));
set(ax(1),'XTickLabel',[]);
set(ax(2),'XTickLabel',[],'YAxisLocation','Right');
ylabel(ax(1),'RL_{min} cm');
ylabel(ax(2),'r_{1,max} cm');
ylabel(ax(3),'R_{1,2}');
grid(ax(3),'on');
xlabel(ax(3),'Phi radians');
linkaxes(ax, 'x');
% xlim(ax(3),minmax(k_set));
%%
figure;
% v = Summary.R1 < 200 & Summary.R2 > -200;
R1 = [Summary.R1];
R2 = [Summary.R2];
colorfield = 'Phi';
m = [Summary.(colorfield)];
v = R1 > 0;
scatter(R1(v),R2(v),[],m(v));
xl = xlim;
R1max = xl(2);
hold on;
h1 = plot([0, R1max]', [L, L-R1max; L L]','k');
yl = ylim;
h2 = plot([L L],yl,'k');
set([h1' h2],'Color',.85*[1 1 1]);
title(sprintf('L = %d', L));
xlabel('R1 cm'); ylabel('R2 cm');
h = colorbar;
title(h,colorfield);
