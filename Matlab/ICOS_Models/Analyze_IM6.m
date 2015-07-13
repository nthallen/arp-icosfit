%%
P = ICOS_Model6.props;
P.visible = 0;
P.HR = 0;
P.visibility = [];
P.focus = 1;
P.R2 = -28.12;
P.y0 = 3.0;
P.z0 = -0.15;
P.dy = 0.097817789473684; % circular for -28.12 which is
P.dz = 0.005696000000000;
P.detector_spacing = 5.105263157894736;
P.Lenses = { 'ZC_PM_25_100' };
P.evaluate_endpoints = -1;
P.injection_scale = 2/3; %1.0;
IB = ICOS_beam(@ICOS_Model6, P);
%%
IB.Sample('opt_n', 5, 'ICOS_passes', 10000);
IB.Integrate;
save IB6_150707_10000x100.mat IB
%%
%IB.Animate;
%%
%%
% P.ICOS_passes_per_injection = 10000; % 10000 for final measurement
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
% opt_n = 5;
% TStart = tic;
% for i = 1:Nsamples
%   if ResLen > NRes
%     P.beam_dy = Dy(i);
%     P.beam_dz = Dz(i);
%     PM = ICOS_Model6(P);
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
% P.beam_dy = 0;
% P.beam_dz = 0;
% save AIM6_150704_10000x100.mat ResE ResD ResP ResNPass ResSample NRes opt_n NPI Nsamples PM
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
%     v = ResNPass == j;
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
%   %Dyz = .01; % detector width, cm
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
%       Pimg(yi,zi) = sum(ResP(Ei(vy & vz)))/(Pexp*Nsamples);
%     end
%   end
%   %%
%   figure;
%   h = image(Yr*dyz,Zr*dyz,Pimg','CDataMapping','scaled'); shg
%   set(gca,'DataAspectRatio',[1 1 1],'Ydir','normal');
%   colorbar;
%   title(sprintf('Fraction of expected total power, %.1fmm Detector', Dyz*10));
%   xlabel('Y position of detector center, cm');
%   ylabel('Z position of detector center, cm');
% end
