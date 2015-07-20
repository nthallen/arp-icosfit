classdef ICOS_beam < handle
  % ICOS_beam provides another layer of modeling
  % to look at the effects of the finite beam
  % dimensions. This class should work with any
  % class derived from opt_model_p. The properties
  % for the model must include:
  %   y0, z0, dy, dz, beam_dy, beam_dz, beam_diameter
  %   T (trasmittance)
  properties
    model % function pointer to the opt_model_p
    P % base properties
    IBP % ICOS_beam properties
    Res % Results structure
  end
  
  methods
    function IB = ICOS_beam(om, P)
      IB.model = om;
      IB.P = P;
      IB.IBP.ICOS_passes = 100;
      IB.IBP.NPI = IB.IBP.ICOS_passes/2;
      IB.IBP.beam_samples = 100;
      IB.IBP.dyz = 0.01;
      IB.IBP.Dyz = [0.01 0.2];
      IB.IBP.opt_n = 0;
      IB.IBP.Tinterval = 60;
      IB.IBP.Track_Power = 0;
      IB.Res = [];
    end
    
    function Sample(IB, varargin)
      if isstruct(IB.Res)
        error('MATLAB:HUARP:Sampled', 'Sample has already been run');
      end
      i = 1;
      while i < length(varargin)
        if isfield(IB.IBP, varargin{i})
          IB.IBP.(varargin{i}) = varargin{i+1};
        elseif isfield(IB.P, varargin{i})
          IB.P.(varargin{i}) = varargin{i+1};
        else
          error('MATLAB:HUARP:InvalidOption', ...
            'Invalid option: "%s"', varargin{i});
        end
        i = i+2;
      end
      IB.P.ICOS_passes_per_injection = IB.IBP.ICOS_passes;
      IB.IBP.NPI = IB.P.ICOS_passes_per_injection/2;
      IB.P.max_rays = floor(IB.P.ICOS_passes_per_injection*5);
      Nsamples = IB.IBP.beam_samples;
      X = rand(Nsamples,1);
      r = IB.P.beam_diameter*sqrt(-log(X))/2;
      th = 2*pi*rand(Nsamples,1);
      Dy = r.*cos(th);
      Dz = r.*sin(th);
      ResLen = Nsamples*IB.IBP.NPI;
      IB.Res.N = 0;
      IB.Res.D = zeros(ResLen, 3);
      IB.Res.E = zeros(ResLen, 3);
      IB.Res.P = zeros(ResLen,1);
      IB.Res.NPass = zeros(ResLen,1);
      IB.Res.Sample = zeros(ResLen,1);
      IB.Res.NPasses = zeros(Nsamples,1);
      IB.Res.perimeter = [];
      IB.Res.Int(length(IB.IBP.Dyz)).img = [];
      IB.Res.Dy = Dy;
      IB.Res.Dz = Dz;
      if IB.IBP.Track_Power
        IB.Res.Pwr = [];
      end
      P = IB.P;
      % IB.IBP.evaluate_endpoints = 0;
      TStart = tic;
      Treport = 0;
      for i = 1:Nsamples
        if ResLen > IB.Res.N
          P.beam_dy = Dy(i);
          P.beam_dz = Dz(i);
          PM = IB.model(P);
          if IB.IBP.opt_n < 1 || IB.IBP.opt_n > length(PM.M.Optic)
            error('MATLAB:HUARP:OpticOOR', 'Invalid Optic Number');
          end
          if isempty(IB.Res.perimeter)
            IB.Res.perimeter = PM.M.Optic{IB.IBP.opt_n}.Surface{1}.perimeter;
          end
          Ri = 2:PM.M.n_rays;
          pre_opt = [PM.M.Rays([PM.M.Rays(Ri).n_inc]).n_opt];
          cur_opt = [PM.M.Rays(Ri).n_opt];
          vo = pre_opt == 2 & cur_opt == 3;
          pass = cumsum(vo);
          % but reset the pass count 
          vr = pre_opt == 1 & cur_opt == 2;
          repass = vr .* pass;
          vri = find(vr);
          repass(vri(2:end)) = diff(repass(vri));
          pass = cumsum(vo-repass)';
          
          vf = find([PM.M.Rays(Ri).n_opt] == IB.IBP.opt_n)+1;
          nvf = length(vf);
          IB.Res.NPasses(i) = nvf;
          if IB.Res.N+nvf > ResLen
            nvf = ResLen - IB.Res.N;
          end
          if nvf > IB.IBP.NPI
            nvf = IB.IBP.NPI;
          end
          for j=1:nvf
            IB.Res.N = IB.Res.N+1;
            ray = PM.M.Rays(vf(j)).ray;
            IB.Res.E(IB.Res.N,:) = ray.E;
            IB.Res.D(IB.Res.N,:) = ray.D;
            IB.Res.P(IB.Res.N) = ray.P;
            IB.Res.NPass(IB.Res.N) = pass(vf(j));
            IB.Res.Sample(IB.Res.N) = i;
          end
          if IB.IBP.Track_Power
            A = [pre_opt' cur_opt'];
            [B,~,ib] = unique(A,'rows');
            ib(pass > IB.IBP.NPI) = 0;
            RP = zeros(size(Ri'));
            inside = zeros(size(Ri'));
            for Pi=1:length(Ri)
              ray = PM.M.Rays(Ri(Pi)).ray;
              RP(Pi) = ray.P;
              inside(Pi) = ray.Inside;
            end
            outside = ~inside;
            for Pi=1:length(B)
              fld = sprintf('R%d_%d', B(Pi,1), B(Pi,2));
              if ~isfield(IB.Res.Pwr, fld)
                IB.Res.Pwr.(fld) = struct('I',0,'O',0,'NI',0,'NO',0);
              end
              for pre = 'IO'
                if pre == 'I'
                  icond = inside;
                else
                  icond = outside;
                end
                IB.Res.Pwr.(fld).(pre) = IB.Res.Pwr.(fld).(pre) + ...
                  sum(RP((ib == Pi) & icond));
                IB.Res.Pwr.(fld).(['N' pre]) = IB.Res.Pwr.(fld).(['N' pre]) + ...
                  sum((ib == Pi) & icond);
              end
            end
            if isfield(IB.Res.Pwr,'R2_3') && isfield(IB.Res.Pwr,'R3_4') ...
                && IB.Res.Pwr.R2_3.NI ~= IB.Res.Pwr.R3_4.NI
              fprintf(1, 'Sample %d: R2_3.NI ~= R3_4.NI\n', i);
            end
          end
        end
        TIter = toc(TStart);
        if TIter > Treport + IB.IBP.Tinterval
          fprintf(1,'%.1f: Iteration: %d Passes: %d Total: %d\n', TIter, ...
            i, IB.Res.NPasses(i), IB.Res.N);
          Treport = TIter;
        end
      end
      T = P.T;
      IB.Res.Ptotal = sum(IB.Res.P(1:IB.Res.N))/Nsamples;
      Pin = T;
      fprintf(1,'Total power is %.1f%% of injected power\n', ...
        100*IB.Res.Ptotal/Pin);
      IB.Res.Pexp = Pin/2 * (1-(1-T)^(IB.IBP.NPI*2));
      fprintf(1,'Total power is %.1f%% of expected power\n', ...
        100*IB.Res.Ptotal/IB.Res.Pexp);
      fprintf(1,'Total power:    %.4g\n', IB.Res.Ptotal);
      fprintf(1,'Expected power: %.4g\n', IB.Res.Pexp);
    end
    
    function Animate(IB, runafter)
      if ~isstruct(IB.Res)
        error('MATLAB:HUARP:Unsampled', ...
          'Must run Sample before Animate()');
      end
      if nargin < 2
        runafter = 50;
      end
      figure;
      NPI = min(50, IB.IBP.NPI);
      Ei = 1:IB.Res.N;
      h = [];
      p = IB.Res.perimeter;
      Yr = minmax([IB.Res.E(Ei,2);p(:,2)]');
      Zr = minmax([IB.Res.E(Ei,3);p(:,3)]');
      plot(p(:,2),p(:,3),'k');
      hold on;
      for j=1:NPI; % NPI
        v = IB.Res.NPass(Ei) == j;
        if ~isempty(h)
          set(h,'MarkerEdgeColor',[0 1 0]);
        end
        h = plot(IB.Res.E(Ei(v),2), IB.Res.E(Ei(v),3), '.r');
        set(gca,'xlim',Yr,'ylim',Zr,'DataAspectRatio',[1 1 1]);
        title(sprintf('%d Passes', j));
        if j > runafter
          drawnow;
          shg;
          pause(0.1);
        else
          pause;
        end
      end
      hold off;
    end
    
    function ff = Integrate(IB)
      if ~isstruct(IB.Res)
        error('MATLAB:HUARP:Unsampled', ...
          'Must run Sample before Integrate()');
      end
      Nsamples = IB.IBP.beam_samples;
      Tstart = tic;
      Treport = 0;
      dyz = IB.IBP.dyz; % spacing of detector centers
      figs = zeros(length(IB.IBP.Dyz),1);
      for Dyzi = 1:length(IB.IBP.Dyz)
        Dyz = IB.IBP.Dyz(Dyzi);
        if isempty(IB.Res.Int(Dyzi).img)
          %Dyz = .01; % detector width, cm
          Ei = 1:IB.Res.N;
          Yr = minmax(IB.Res.E(Ei,2)');
          Yr(1) = floor((Yr(1)-Dyz/2)/dyz);
          Yr(2) = ceil((Yr(2)+Dyz/2)/dyz);
          Zr = minmax(IB.Res.E(Ei,2)');
          Zr(1) = floor((Zr(1)-Dyz/2)/dyz);
          Zr(2) = ceil((Zr(2)+Dyz/2)/dyz);
          IB.Res.Int(Dyzi).Yr = Yr;
          IB.Res.Int(Dyzi).Zr = Zr;
          Y = dyz*(Yr(1):Yr(2));
          Z = dyz*(Zr(1):Zr(2));
          nY = length(Y);
          nZ = length(Z);
          Pimg = zeros(nY,nZ);
          Ptotal = sum(IB.Res.P(Ei))/Nsamples;
          for yi = 1:nY
            TIter = toc(Tstart);
            if TIter - Treport > 10
              fprintf(1, 'T = %.1f Dyz = %.2f yi = %d/%d\n', TIter, Dyz, yi, nY);
              Treport = TIter;
            end
            dy = IB.Res.E(Ei,2)-Y(yi);
            vy = dy > -Dyz/2 & dy <= Dyz/2;
            for zi = 1:nZ
              dz = IB.Res.E(Ei,3)-Z(zi);
              vz = dz > -Dyz/2 & dz <= Dyz/2;
              Pimg(yi,zi) = sum(IB.Res.P(Ei(vy & vz)))/(IB.Res.Pexp*Nsamples);
            end
          end
          IB.Res.Int(Dyzi).img = Pimg;
        else
          Pimg = IB.Res.Int(Dyzi).img;
          Yr = IB.Res.Int(Dyzi).Yr;
          Zr = IB.Res.Int(Dyzi).Zr;
        end
        figs(Dyzi) = figure;
        h = image(Yr*dyz,Zr*dyz,Pimg','CDataMapping','scaled'); shg
        pos = get(figs(Dyzi),'position');
        pos(1) = 5 + (Dyzi-1)*pos(3);
        set(figs(Dyzi),'position',pos);
        set(gca,'DataAspectRatio',[1 1 1],'Ydir','normal');
        ch = colorbar;
        mname = strrep(func2str(IB.model),'_','\_');
        
        title(sprintf('%s: %.1fmm Detector %d passes %d samples', ...
          mname, Dyz*10, IB.IBP.ICOS_passes, ...
          IB.IBP.beam_samples ));
        ylabel(ch,'Fraction of expected total power');
        xlabel('Y position of detector center, cm');
        ylabel('Z position of detector center, cm');
        drawnow; shg;
      end
      if nargout > 0
        ff = figs;
      end
    end
    
    function display(IB)
      fprintf(1,'ICOS_beam:\n');
      fprintf(1,'  model: %s\n', func2str(IB.model));
      fprintf(1,'  ICOS_passes: %d\n', IB.IBP.ICOS_passes);
      fprintf(1,'  beam_samples: %d\n', IB.IBP.beam_samples);
      fprintf(1,'  dyz: %.4f\n', IB.IBP.dyz);
      fprintf(1,'  Dyz: [');
      fprintf(1,' %.3f', IB.IBP.Dyz);
      fprintf(1,' ]\n');
      if isstruct(IB.Res)
        if isempty(IB.Res.Int(1).img)
          fprintf(1,'  Sampled, Unintegrated\n');
        else
          fprintf(1,'  Sampled, Integrated\n');
        end
      else
        fprintf(1,'  Unsampled\n');
      end
    end
    
    function draw(IB, varargin)
      P = IB.P;
      P.visible = 1;
      % P.visibility = [];
      i = 1;
      while i < length(varargin)
        if isfield(P, varargin{i})
          P.(varargin{i}) = varargin{i+1};
        else
          error('MATLAB:HUARP:InvalidOption', ...
            'Invalid option: "%s"', varargin{i});
        end
        i = i+2;
      end
      PM = IB.model(P);
    end
    
    function png(IB, bn)
      % IB.png(base)
      % Run IB.Integrate and saves images as .png files
      % named <base>_<ICOS_passes>x<beam_samples>_<Dyz>.png
      if nargin < 2
        bn = func2str(IB.model);
      end
      figs = IB.Integrate;
      for Dyzi = 1:length(IB.IBP.Dyz)
        fname = sprintf('%s_%dx%d_%.2f.png', bn, IB.IBP.ICOS_passes, ...
          IB.IBP.beam_samples, IB.IBP.Dyz(Dyzi));
        print(figs(Dyzi), '-dpng', fname);
        delete(figs(Dyzi));
        fprintf(1,'Output written to %s\n', fname);
      end
    end
  end
end
