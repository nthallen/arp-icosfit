%%
% This file will perform Monte Carlo simulation of a gaussian
% beam in order to better approximate power on the final detector.
%
% These are the final parameters drawn from Explore_IM5.m
% This describes a cylindrical ICOS configuration with
% Mirrors with RC of 75 cm spaced at 50 cm.
P = ICOS_Model5.props;
P.visible = 0;
P.visibility = [];
P.RC = 75;
P.focus = 1;
P.herriott_spacing = 12.8;
P.plot_endpoints = 0;
P.evaluate_endpoints = -1;
P.fmt_HRC = 'Herriott RoC = %.2f';
P.pause = 0;
P.stop_ICOS = 0;
P.HR = 0;
P.ICOS_passes_per_injection = 63;
P.max_rays = 3000;
P.HRC = 40.64;
P.y0 = 3.0;
P.z0 = 0;
P.injection_scale = 2/3;
P.dy = 0.097696842105263; % circular
P.dz = 0.057006315789474;
P.Lenses = {'Lens1', 'ZC_NM_25_100' };
P.Lens_Space = [0.2, 6.8];
P.detector_spacing = 2.308;
%%
IB = ICOS_beam(@ICOS_Model5,P);
IB.Sample('opt_n', 6, 'ICOS_passes', 100);
IB.Integrate;
%IB.png('IB5r3');
%save IB5r3_150707_10000x100.mat IB
%%
IB.Animate
%% These two sections are for displaying a single pass of all the samples
P.beam_dy = IB.Res.Dy(1);
P.beam_dz = IB.Res.Dz(1);
P.visible = 1;
PM = ICOS_Model5(P);
%%
m = P.injection_scale;
for i = 2:length(IB.Res.Dy)
  % disp(i);
  Pincident = PM.M.Optic{2}.O + [0, m*P.y0 + IB.Res.Dy(i), m*P.z0 + IB.Res.Dz(i)];
  Dincident = [1, -m*P.dy, m*P.dz];
  Oincident = Pincident - (P.herriott_spacing+1)*Dincident;
  PM.M.push_ray(opt_ray(Oincident, Dincident), 0, 1, 2);
  PM.M.propagate;
end
%%
%PM = ICOS_Model5(P);
%%
%P.visibility = [0 0 0 0 0];
%PM = ICOS_Model5(P);

%% Now lets see what we can do with monte carlo sampling
% First, let's back back down to just ICOS
% How many passes should we look at?
%   We need to collect light from about 1000 passes in order to count
%   about 90%. To reach 99% takes about 18000 bounces.
% P.HR = 0;
% P.visibility = [];
% P.visibility = [0 0 0 0 0];
% P.visible = 0;
% P.draw_endpoints = 0;
% P.evaluate_endpoints = 0;
% P.ICOS_passes_per_injection = 1000;
% P.max_rays = floor(P.ICOS_passes_per_injection*5);
% Nsamples = 100;
% NPI = P.ICOS_passes_per_injection/2;
% X = rand(Nsamples,1);
% r = P.beam_diameter*sqrt(-log(X))/2;
% th = 2*pi*rand(Nsamples,1);
% Dy = r.*cos(th);
% Dz = r.*sin(th);
% ResLen = Nsamples*NPI;
% NRes = 0;
% ResD = zeros(ResLen, 3);
% ResE = ResD;
% ResP = zeros(ResLen,1);
% ResNPass = zeros(ResLen,1);
% ResSample = zeros(ResLen,1);
% NPasses = zeros(Nsamples,1);
% y0 = P.y0;
% z0 = P.z0;
% opt_n = 6;
% TStart = tic;
% for i = 1:Nsamples
%   if ResLen > NRes
%     P.y0 = y0+Dy(i);
%     P.z0 = z0+Dz(i);
%     PM = ICOS_Model5(P);
%     if opt_n < 1 || opt_n > length(PM.M.Optic)
%       error('MATLAB:HUARP:OpticOOR', 'Invalid Optic Number');
%     end
%     vf = find([PM.M.Rays(1:PM.M.n_rays).n_opt] == opt_n);
%     nvf = length(vf);
%     NPasses(i) = nvf;
%     if NRes+nvf > ResLen
%       nvf = ResLen - NRes;
%     end
%     if nvf > NPI
%       nvf = NPI;
%     end
%     for j=1:nvf
%       NRes = NRes+1;
%       ResE(NRes,:) = PM.M.Rays(vf(j)).ray.E;
%       ResD(NRes,:) = PM.M.Rays(vf(j)).ray.D;
%       ResP(NRes) = PM.M.Rays(vf(j)).ray.P;
%       ResNPass(NRes) = j;
%       ResSample(NRes) = i;
%     end
%   end
%   TIter = toc(TStart);
%   fprintf(1,'%.1f: Iteration: %d NRes: %d\n', TIter, i, NRes);
% end
% if NRes ~= ResLen
%   fprintf(1,'NRes = %d, ResLen = %d\n', NRes, ResLen);
% end
% P.y0 = y0;
% P.z0 = z0;
% save AIM5_150705_1000x100.mat ResE ResD ResP ResNPass ResSample NRes opt_n NPI Nsamples PM
% %%
% T = 250e-6;
% Ptotal = sum(ResP(1:NRes))/Nsamples;
% Pin = T;
% fprintf(1,'Total power is %.1f%% of injected power\n', 100*Ptotal/Pin);
% Pexp = Pin/2 * (1-(1-T)^(NPI*2));
% fprintf(1,'Total power is %.1f%% of expected power\n', 100*Ptotal/Pexp);
% fprintf(1,'Total power:    %.4g\n', Ptotal);
% fprintf(1,'Expected power: %.4g\n', Pexp);
% %%
% if NPI < 100
%   %%
%   figure;
%   Ei = 1:NRes;
%   h = [];
%   p = PM.M.Optic{opt_n}.Surface{1}.perimeter;
%   Yr = minmax([ResE(Ei,2);p(:,2)]');
%   Zr = minmax([ResE(Ei,3);p(:,3)]');
%   plot(p(:,2),p(:,3),'k');
%   hold on;
%   for j=1:50; % NPI
%     v = mod(Ei-1,NPI) == j-1;
%     if ~isempty(h)
%       set(h,'MarkerEdgeColor',[0 1 0]);
%     end
%     h = plot(ResE(Ei(v),2), ResE(Ei(v),3), '.r');
%     set(gca,'xlim',Yr,'ylim',Zr,'DataAspectRatio',[1 1 1]);
%     title(sprintf('%d Passes', j));
%     pause;
%   end
%   hold off;
% end
% %%
% for Dyz = [.01 .2]
%   dyz = .01; % spacing of detector centers
%   % Dyz = .2; % detector width, cm
%   Ei = 1:NRes;
%   Yr = minmax(ResE(Ei,2)');
%   Yr(1) = floor((Yr(1)-Dyz/2)/dyz);
%   Yr(2) = ceil((Yr(2)+Dyz/2)/dyz);
%   Zr = minmax(ResE(Ei,2)');
%   Zr(1) = floor((Zr(1)-Dyz/2)/dyz);
%   Zr(2) = ceil((Zr(2)+Dyz/2)/dyz);
%   Y = dyz*(Yr(1):Yr(2));
%   Z = dyz*(Zr(1):Zr(2));
%   nY = length(Y);
%   nZ = length(Z);
%   Pimg = zeros(nY,nZ);
%   Ptotal = sum(ResP(1:NRes))/Nsamples;
%   for yi = 1:nY
%     fprintf(1, 'Dyz = %.2f yi = %d/%d\n', Dyz, yi, nY);
%     dy = ResE(Ei,2)-Y(yi);
%     vy = dy > -Dyz/2 & dy <= Dyz/2;
%     for zi = 1:nZ
%       dz = ResE(Ei,3)-Z(zi);
%       vz = dz > -Dyz/2 & dz <= Dyz/2;
%       Pimg(yi,zi) = sum(ResP(Ei(vy & vz)))/(Ptotal*Nsamples);
%     end
%   end
%   %%
%   figure;
%   h = image(Yr*dyz,Zr*dyz,Pimg','CDataMapping','scaled'); shg
%   set(gca,'DataAspectRatio',[1 1 1],'Ydir','normal');
%   colorbar;
%   title(sprintf('Fraction of total expected power, %.1fmm Detector', Dyz*10));
%   xlabel('Y position of detector center, cm');
%   ylabel('Z position of detector center, cm');
% end
%% How many different rays to be representative? I'm inclined to say 1000
% N = 1000;
% X = rand(100,1);
% w = 0.2;
% r = w*sqrt(-log(X));
% [n,x] = hist(r,20);
% dx = mean(diff(x));
% fn = n./(N*dx*2*pi*x);
% plot(x,fn,'o'); shg;
% %%
% plot(x,n); shg