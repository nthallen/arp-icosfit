classdef ICOS_beam < handle
  % ICOS_beam provides another layer of modeling
  % to look at the effects of the finite beam
  % dimensions. This class should work with any
  % class derived from opt_model_p. The properties
  % for the model must include:
  %   y0, z0, dy, dz, beam_dy, beam_dz, beam_diameter
  %   T (trasmittance)
  % See also: ICOS_beam.ICOS_beam, ICOS_beam.Sample,
  % ICOS_beam.PowerSummary, ICOS_beam.Animate, ICOS_beam.Integrate,
  % ICOS_beam.display, ICOS_beam.draw, ICOS_beam.png
  properties
    model % function pointer to the opt_model_p
    P % base properties
    IBP % ICOS_beam properties
    Res % Results structure
  end
  
  methods
    function IB = ICOS_beam(om, P, varargin)
      % IB = ICOS_beam(model, P)
      % IB.IBP contains ICOS_beam options, including:
      %   Herriott_passes: advisory on the maximum number
      %    of passes expected in the Herriott cell. Used to
      %    reserve space for the rays
      %   n_optics: advisory on the number of optics in the model.
      %    This will be tested after the model is run.
      IB.model = om;
      IB.P = P;
      IB.IBP.mnc = 'ib.def';
      IB.IBP.Herriott_passes = 100;
      IB.IBP.ICOS_passes = 100;
      IB.IBP.theta = 15; % acceptance angle
      IB.IBP.n_optics = 5;
      IB.IBP.NPI = IB.IBP.ICOS_passes;
      IB.IBP.beam_samples = 100;
      IB.IBP.beam_divergence = 0;
      IB.IBP.dyz = 0.01;
      IB.IBP.Dyz = [0.01 0.2];
      IB.IBP.opt_n = 0;
      IB.IBP.Tinterval = 60;
      IB.IBP.Int.Yr = []; % 2x2, row for each Dyz
      IB.IBP.Int.Zr = [];
      IB.IBP.Track_Power = 0;
      IB.IBP.Npass_dim = 2; % The number of different pass counts (Herriott & ICOS)
      IB.IBP.rng_state = [];
      IB.Res = [];
      % IB.Res includes information about all rays incident on optic opt_n
      % IB.Res.N: total number of rays recorded
      % IB.Res.D: Nx3 direction vectors
      % IB.Res.E: Nx3 endpoint vectors
      % IB.Res.P: Nx1 power
      % IB.Res.NPass: N x Npass_dim: which pass
      % IB.Res.Sample: N x 1: which sample
      % IB.Res.NPasses: beam_samples x 1: number of rays incident on optic opt_n for each sample
      % IB.Res.perimeter: np x 3: outline of optic opt_n
      % IB.Res.Int(length(Dyz)).img: image of power integration
      % IB.Res.Dy, Dz: beam_samples: relative position of beam sample
      % IB.Res.Pwr: Summary of power for each transition
      % IB.Res.Ptotal: Total power (by some metric)
      % IB.Res.Pexp: Expected power (by some possibly different metric)
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
    end
    
    function Sample(IB, varargin)
      % IB.Sample(options)
      % options are pairs of arguments as keyword, value. These may include
      % any fields in the IB.IBP structure or any fields in the model
      % parameter structure IB.P.
      %
      % IB.IBP options:
      %   opt_n: The index of the optic whose rays we want to track
      %   Track_Power: boolean: default is 0
      %   ICOS_passes: The number of ICOS passes to track. Default is 100
      %   beam_samples: The number of random beams. Default is 100
      %
      % IB.Sample generates a random sample of rays to represent the
      % physical dimensions of a real laser beam and propagates the optical
      % model for each, saving the results within the IB object.
      % See also: ICOS_beam, ICOS_beam.Integrate
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
      if ~isempty(IB.IBP.rng_state)
        rng(IB.IBP.rng_state);
      else
        IB.IBP.rng_state = rng;
      end
      IB.P.ICOS_passes_per_injection = IB.IBP.ICOS_passes;
      IB.IBP.NPI = IB.P.ICOS_passes_per_injection;
      IB.P.max_rays = floor(IB.P.ICOS_passes_per_injection* ...
        IB.IBP.n_optics * IB.IBP.Herriott_passes);
      Nsamples = IB.IBP.beam_samples;
      X = rand(Nsamples,1);
      r = IB.P.beam_diameter*sqrt(-log(X))/2;
      th = 2*pi*rand(Nsamples,1);
      Dy = r.*cos(th);
      Dz = r.*sin(th);
      ResLen = Nsamples*IB.IBP.NPI*IB.IBP.Herriott_passes;
      IB.Res.N = 0;
      IB.Res.D = zeros(ResLen, 3);
      IB.Res.E = zeros(ResLen, 3);
      IB.Res.P = zeros(ResLen,1);
      IB.Res.NPass = zeros(ResLen,IB.IBP.Npass_dim);
      IB.Res.Sample = zeros(ResLen,1);
      IB.Res.NPasses = zeros(Nsamples,1);
      IB.Res.perimeter = [];
      IB.Res.Int(length(IB.IBP.Dyz)).img = [];
      IB.Res.Dy = Dy;
      IB.Res.Dz = Dz;
      IB.Res.Warn.n_rays = 0;
      IB.Res.Warn.resovf = 0;
      IB.Res.Warn.resovf2 = 0;
      if IB.IBP.Track_Power
        IB.Res.Pwr = [];
      end
      P = IB.P;
      % IB.IBP.evaluate_endpoints = 0;
      TStart = tic;
      Treport = 0;
      for i = 1:Nsamples
        if ResLen <= IB.Res.N
          warning('Matlab:HUARP:SampleSkip', ...
            'Skipping at sample %d of %d: result array inadequate', ...
            i, Nsamples);
          break;
        else
          P.beam_dy = Dy(i);
          P.beam_dz = Dz(i);
          PM = IB.model(P);
          if IB.IBP.opt_n < 1 || IB.IBP.opt_n > length(PM.M.Optic)
            error('MATLAB:HUARP:OpticOOR', 'Invalid Optic Number');
          end
          if IB.IBP.n_optics ~= length(PM.M.Optic)
            warning('Matlab:HUARP:n_optics', ...
              'IB.IBP.n_optics (%d) does not match model (%d)', ...
              IB.IBP.n_optics, length(PM.M.Optic));
          end
          if PM.M.n_rays >= PM.M.max_rays
            warning('Matlab:HUARP:model_maxed', ...
              'n_rays >= max_rays (maxed out): accuracy will be limited');
            IB.Res.Warn.n_rays = IB.Res.Warn.n_rays+1;
          end
          if isempty(IB.Res.perimeter)
            IB.Res.perimeter = PM.M.Optic{IB.IBP.opt_n}.Surface{1}.perimeter;
          end
          Ri = 1:PM.M.n_rays;
          
          vf = find([PM.M.Rays(Ri).n_opt] == IB.IBP.opt_n);
          nvf = length(vf);
          IB.Res.NPasses(i) = nvf;
          if IB.Res.N+nvf > ResLen
            warning('Matlab:HUARP:ResultOverflow', ...
              'Result overflow in sample %d of %d', i, Nsamples);
            nvf = ResLen - IB.Res.N;
            IB.Res.Warn.resovf = IB.Res.Warn.resovf+1;
          end
          if nvf > IB.IBP.NPI*IB.IBP.Herriott_passes
            warning('Matlab:HUARP:ResultOverflow', ...
              'nvf > NPI*HP in sample %d of %d', i, Nsamples);
            IB.Res.Warn.resovf2 = IB.Res.Warn.resovf2+1;
            % nvf = IB.IBP.NPI*IB.IBP.Herriott_passes;
          end
          % Save rays associated with IB.IBP.opt_n
          for j=1:nvf
            IB.Res.N = IB.Res.N+1;
            ray = PM.M.Rays(vf(j)).ray;
            IB.Res.E(IB.Res.N,:) = ray.E;
            IB.Res.D(IB.Res.N,:) = ray.D;
            IB.Res.P(IB.Res.N) = ray.P;
            NPd = min(length(ray.pass),IB.IBP.Npass_dim);
            IB.Res.NPass(IB.Res.N,1:NPd) = ray.pass(1:NPd);
            IB.Res.Sample(IB.Res.N) = i;
          end
          if IB.IBP.Track_Power
            % Summarize power associated with all optics
            pre_opt = [0,PM.M.Rays([PM.M.Rays(Ri(2:end)).n_inc]).n_opt];
            cur_opt = [PM.M.Rays(Ri).n_opt];
            A = [pre_opt' cur_opt'];
            [B,~,ib] = unique(A,'rows');
            % ib(pass > IB.IBP.NPI) = 0;
            RP = zeros(size(Ri'));
            inside = zeros(size(Ri'));
            outside = zeros(size(Ri'));
            ap_exit = zeros(size(Ri'));
            % Categorize all the rays
            for Pi=1:length(Ri)
              ray = PM.M.Rays(Ri(Pi)).ray;
              RP(Pi) = ray.P;
              inside(Pi) = ray.Inside > 0;
              ap_exit(Pi) = ray.Inside < 0;
              outside(Pi) = ray.Inside == 0;
            end
            % Now loop through the transitions
            for Pi=1:size(B,1)
              fld = sprintf('R%d_%d', B(Pi,1), B(Pi,2));
              if ~isfield(IB.Res.Pwr, fld)
                if strcmp(fld, 'R2_1')
                  IB.Res.Pwr.(fld) = struct('I',0,'O',0,'NI',0,'NO',0, ...
                    'E1',0,'NE1',0,'EO',0,'NEO',0);
                else
                  IB.Res.Pwr.(fld) = struct('I',0,'O',0,'NI',0,'NO',0);
                end
              end
              transitionV = ib == Pi; % boolean of rays that match trans.
              for pre = 'IO'
                if pre == 'I'
                  icond = inside;
                else
                  icond = outside;
                end
                IB.Res.Pwr.(fld).(pre) = IB.Res.Pwr.(fld).(pre) + ...
                  sum(RP(transitionV & icond));
                IB.Res.Pwr.(fld).(['N' pre]) = IB.Res.Pwr.(fld).(['N' pre]) + ...
                  sum(transitionV & icond);
              end
              if strcmp(fld, 'R2_1')
                ei = find(ap_exit,1);
                if ~isempty(ei)
                  if sum(transitionV) > 1
                    IB.Res.Pwr.R2_1.EO = IB.Res.Pwr.R2_1.EO + RP(ei);
                    IB.Res.Pwr.R2_1.NEO = IB.Res.Pwr.R2_1.NEO+1;
                  else
                    IB.Res.Pwr.R2_1.E1 = IB.Res.Pwr.R2_1.E1 + RP(ei);
                    IB.Res.Pwr.R2_1.NE1 = IB.Res.Pwr.R2_1.NE1+1;
                  end
                end
              end
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
    end
    
    function Pwr = PowerSummary(IB)
      % IB.PowerSummary;
      % Pwr = IB.PowerSummary;
      %
      % Reports on how optical power propagates through the setup.
      % If an output argument is provided, the verbose output is
      % suppressed.
      %
      % The cumulative power loss figure is expected to be a better
      % predictor of performance than the total power ratios.
      %
      % Total power as a percent of injected power is mostly an indication
      % of how many ICOS passes were simulated..
      %
      % Total power as a percent of expected power normalizes for the
      % number of ICOS passes, but yields a number somewhat higher
      % than what we expect as the number of ICOS passes increases.
      % This is because once the beams leave the ICOS cell, they no
      % longer propagate.
      %
      % The Herriott cell details do not add up to the total Herriott cell
      % loss, but should be useful as a general guide as to how the
      % reinjection module is working. The actual reduction in effective
      % power when a beam leaves the Herriott cell depends on how many
      % passes it has already traversed. The detail simply reports the
      % amount of power leaving the cell by various means, but a beam
      % exiting on pass 32 of 33 represents a smaller net loss than one
      % exiting on the first bounce.
      %
      % See also: ICOS_beam, ICOS_beam.display, ICOS_beam.draw,
      % ICOS_beam.png
      if ~isstruct(IB.Res)
        error('MATLAB:HUARP:Unsampled', ...
          'Must run Sample before PowerSummary()');
      end
      T = IB.P.T;
      Nsamples = IB.IBP.beam_samples;
      IB.Res.Ptotal = sum(IB.Res.P(1:IB.Res.N))/Nsamples;
      % Herriott optimal number of passes
      dNP = [ diff(IB.Res.NPass(1:IB.Res.N,:)); -ones(1,size(IB.Res.NPass,2))];
      v = find(IB.Res.NPass(1:IB.Res.N,1)>1 & ...
        (dNP(:,1)<0 | (dNP(:,1)==0 & dNP(:,2)<0)));
      if isempty(v)
        NH = 1;
      else
        NH = mode(IB.Res.NPass(v,1));
      end
      RH = IB.P.HR;
      P2_Ideal = (1-RH^NH)/(1-RH); % normalized by Nsamples
      Pin = T*P2_Ideal;
      Pexp = Pin/2 * (1-(1-T)^(IB.IBP.NPI*2));
      if nargout > 0
        Pwr.NH = NH;
      else
        fprintf(1,'%d Herriott passes\n', NH);
        fprintf(1,'%d ICOS passes\n', IB.IBP.NPI);
        fprintf(1,'Total power is %.1f%% of injected power\n', ...
          100*IB.Res.Ptotal/Pin);
        fprintf(1,'Total power is %.1f%% of expected power\n', ...
          100*IB.Res.Ptotal/Pexp);
      end
      if isfield(IB.Res,'Pwr')
        % Herriott Cell Analysis
        if isfield(IB.Res.Pwr,'R1_2')
          P2_Actual = (IB.Res.Pwr.R0_2.I + IB.Res.Pwr.R1_2.I)/Nsamples;
          Hgain = Nsamples*P2_Actual/(IB.Res.Pwr.R0_2.I+IB.Res.Pwr.R0_2.O);
        else
          P2_Actual = IB.Res.Pwr.R0_2.I/Nsamples;
        end
        P2_loss_pct = 100*(P2_Ideal-P2_Actual)/P2_Ideal;
        cum_pwr = 1-P2_loss_pct/100;
        if nargout == 0
          fprintf(1,'Herriott cell loss: %.1f%%\n', P2_loss_pct);
        end
        if isfield(IB.Res.Pwr,'R2_1') && isfield(IB.Res.Pwr,'R1_2') && ...
            isfield(IB.Res.Pwr.R2_1,'NEO')
          PO_E1 = 100*IB.Res.Pwr.R2_1.E1/Nsamples;
          PO_EO = 100*IB.Res.Pwr.R2_1.EO/Nsamples;
          PO_1O = 100*IB.Res.Pwr.R2_1.O/Nsamples;
          PO_2O = 100*(IB.Res.Pwr.R0_2.O + IB.Res.Pwr.R1_2.O)/Nsamples;
          PO_ICOS = 100*(IB.Res.Pwr.R0_2.I + IB.Res.Pwr.R1_2.I)*IB.P.T/Nsamples;
          PO_Hloss = 100*IB.Res.Pwr.R2_1.I * (1-IB.P.HR)/Nsamples;
          if nargout == 0
            fprintf(1, '  %4.1f%% loss through aperature on 1st pass\n', PO_E1);
            fprintf(1, '  %4.1f%% loss outside on Herriott mirror\n', PO_1O);
            fprintf(1, '  %4.1f%% loss outside on back of ICOS mirror\n', PO_2O);
            fprintf(1, '  %4.1f%% reflective loss on Herriott mirror\n', PO_Hloss);
            fprintf(1, '  %4.1f%% power through aperature later\n', PO_EO);
            fprintf(1, '  %4.1f%% power into ICOS cell\n', PO_ICOS);
          else
            Pwr.H.PO_E1 = PO_E1;
            Pwr.H.PO_EO = PO_EO;
            Pwr.H.PO_1O = PO_1O;
            Pwr.H.PO_2O = PO_2O;
            Pwr.H.PO_ICOS = PO_ICOS;
            Pwr.H.PO_Hloss = PO_Hloss;
          end
        end
        % ICOS analysis
        Iin_Actual = P2_Actual * T;
        I_loss_pct = 100*(IB.Res.Pwr.R2_3.O+IB.Res.Pwr.R3_2.O)/ ...
          (Iin_Actual*Nsamples);
        if nargout == 0
          fprintf(1,'ICOS cell loss: %.1f%%\n', I_loss_pct);
        end
        cum_pwr = cum_pwr * (1-I_loss_pct/100);
        n_focus_optics = length(IB.P.Lenses)+1;
        focus_losses = ones(n_focus_optics,1);
        for opt_n = 3+(0:n_focus_optics-1)
          fld = sprintf('R%d_%d', opt_n, opt_n+1);
          if isfield(IB.Res.Pwr, fld)
            abs_loss = IB.Res.Pwr.(fld).O;
            loss_pct = 100*abs_loss/(IB.Res.Pwr.(fld).I+abs_loss);
            focus_losses(opt_n-2) = loss_pct;
            if nargout == 0
              fprintf(1,'Loss at optic %d: %.1f %%\n', opt_n+1, loss_pct);
            end
            det_pwr_summary = IB.Res.Pwr.(fld);
            % cum_pwr = cum_pwr * (1-loss_pct/100);
          else
            warning('MATLAB:HUARP:NoPwr', 'Expected Power field "%s"', ...
              fld);
          end
        end
        for i = 1:n_focus_optics-1
          cum_pwr = cum_pwr * (1-focus_losses(i)/100);
        end
        cum_pwr0 = cum_pwr * (1-focus_losses(end)/100);
        cum_loss_pct0 = 100*(1-cum_pwr0);
        eff_pwr0 = cum_pwr0 * P2_Ideal;
        if nargout == 0
          fprintf(1, 'Cumulative loss with detector centered: %.1f %%\n', cum_loss_pct0);
          fprintf(1, 'Effective output with detector centered: %.1f\n', eff_pwr0);
        end
        if ~isempty(IB.Res.Int(end).img)
          max_pwr = max(max(IB.Res.Int(end).img));
          IntNorm = Nsamples*T*(1-(1-T)^(IB.IBP.NPI*2))/(2-T);
          det_pwr = max_pwr * IntNorm;
          eff = det_pwr/(det_pwr_summary.I+det_pwr_summary.O);
          Dloss = 100*(1-eff);
          cum_pwr1 = cum_pwr * (1-Dloss/100);
          cum_loss_pct1 = 100*(1-cum_pwr1);
          eff_pwr1 = cum_pwr1 * P2_Ideal;
          if nargout == 0
            fprintf(1, 'Output by moving detector: %.1f\n', max_pwr);
            fprintf(1, 'Loss at detector after translation: %.1f %%\n', ...
              Dloss);
            fprintf(1, 'Cumulative loss after translation: %.1f %%\n', cum_loss_pct1);
            fprintf(1, 'Effective output after translation: %.1f\n', eff_pwr1);
          else
            Pwr.cum_loss_pct = cum_loss_pct1;
          end
        else
          Dloss = focus_losses(end);
          eff_pwr1 = eff_pwr0;
          max_pwr = eff_pwr0;
        end
        if nargout > 0
          Pwr.H_loss = P2_loss_pct;
          Pwr.I_loss = I_loss_pct;
          Pwr.D_loss = Dloss;
          Pwr.F_loss = focus_losses(1:end-1);
          Pwr.eff_pwr = eff_pwr1;
          Pwr.max_pwr = max_pwr;
        end
      end
    end
    
    function Animate(IB, varargin)
      % IB.Animate;
      % IB.Animate('runafter', 15);
      %  'IPass', 1
      %  'HPass', 1
      if ~isstruct(IB.Res)
        error('MATLAB:HUARP:Unsampled', ...
          'Must run Sample before Animate()');
      end
      opt.runafter = 50;
      opt.HPass = 1:min(50,max(IB.Res.NPass(:,1)));
      if IB.IBP.Npass_dim > 1
        opt.IPass = 1:min(50,max(IB.Res.NPass(:,2)));
      else
        opt.IPass = 1;
      end
      i = 1;
      while i < length(varargin)
        if isfield(opt, varargin{i})
          opt.(varargin{i}) = varargin{i+1};
        else
          error('MATLAB:HUARP:InvalidOption', ...
            'Animate: Invalid option: "%s"', varargin{i});
        end
        i = i+2;
      end
      figure;
      Ei = 1:IB.Res.N;
      h = [];
      p = IB.Res.perimeter;
      Yr = minmax([IB.Res.E(Ei,2);p(:,2)]');
      Zr = minmax([IB.Res.E(Ei,3);p(:,3)]');
      plot(p(:,2),p(:,3),'k');
      hold on;
      j = 1;
      for Hi = opt.HPass
        for Ii = opt.IPass
          v = IB.Res.NPass(Ei,1) == Hi;
          if IB.IBP.Npass_dim > 1
            v = v & IB.Res.NPass(Ei,2) == Ii;
          end
          if ~isempty(h)
            set(h,'MarkerEdgeColor',[0 1 0]);
          end
          h = plot(IB.Res.E(Ei(v),2), IB.Res.E(Ei(v),3), '.r');
          set(gca,'xlim',Yr,'ylim',Zr,'DataAspectRatio',[1 1 1]);
          title(sprintf('HPass: %d  IPass: %d', Hi, Ii));
          if j > opt.runafter
            drawnow;
            shg;
            pause(0.1);
          else
            pause;
          end
          j = j+1;
        end
      end
      hold off;
    end
    
    function ff = Integrate(IB)
      % ff = IB.Integrate;
      % Analyzes rays from the Sample operation to generate heat
      % maps in the plane of a particular optical element, usually
      % the detector. If an output argument is specified, the
      % heat maps are displayed.
      % All options are set during the Sample phase.
      if ~isstruct(IB.Res)
        error('MATLAB:HUARP:Unsampled', ...
          'Must run Sample before Integrate()');
      end
      if ~isfield(IB.IBP,'mnc')
        IB.IBP.mnc = 'ib.def';
      end
      if ~isfield(IB.IBP,'theta')
        IB.IBP.theta = 0;
      end
      Nsamples = IB.IBP.beam_samples;
      T = IB.P.T;
      % Pexp is power from Nsamples beams propagated
      % with IB.IBP.NPI ICOS passes
      Pexp = Nsamples*T*(1-(1-T)^(IB.IBP.NPI*2))/(2-T);
      Tstart = tic;
      Treport = 0;
      dyz = IB.IBP.dyz; % spacing of detector centers
      if IB.IBP.theta > 0
        figs = zeros(2*length(IB.IBP.Dyz),1);
      else
        figs = zeros(length(IB.IBP.Dyz),1);
      end
      for Dyzi = 1:length(IB.IBP.Dyz)
        Dyz = IB.IBP.Dyz(Dyzi);
        if isempty(IB.Res.Int(Dyzi).img)
          %Dyz = .01; % detector width, cm
          Ei = 1:IB.Res.N;
          if isempty(IB.IBP.Int.Yr)
            Yr = minmax(IB.Res.E(Ei,2)');
          else
            Yr = IB.IBP.Int.Yr(Dyzi,:);
          end
          Yr(1) = floor((Yr(1)-Dyz/2)/dyz);
          Yr(2) = ceil((Yr(2)+Dyz/2)/dyz);
          if isempty(IB.IBP.Int.Zr)
            Zr = minmax(IB.Res.E(Ei,2)');
          else
            Zr = IB.IBP.Int.Zr(Dyzi,:);
          end
          Zr(1) = floor((Zr(1)-Dyz/2)/dyz);
          Zr(2) = ceil((Zr(2)+Dyz/2)/dyz);
          IB.Res.Int(Dyzi).Yr = Yr;
          IB.Res.Int(Dyzi).Zr = Zr;
          Y = dyz*(Yr(1):Yr(2));
          Z = dyz*(Zr(1):Zr(2));
          nY = length(Y);
          nZ = length(Z);
          Pimg = zeros(nY,nZ);
          if IB.IBP.theta > 0
            Rimg = zeros(nY,nZ);
            thOK = atand(sqrt((IB.Res.D(Ei,2)./IB.Res.D(Ei,1)).^2 + ...
              (IB.Res.D(Ei,3)./IB.Res.D(Ei,1)).^2)) <= IB.IBP.theta;
          else
            thOK = ones(size(Ei))>0;
          end
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
              Pimg(yi,zi) = sum(IB.Res.P(Ei(vy & vz & thOK)))/Pexp;
              if IB.IBP.theta > 0
                Rimg(yi,zi) = sum(IB.Res.P(Ei(vy & vz & ~thOK)))/Pexp;
              end
            end
          end
          IB.Res.Int(Dyzi).img = Pimg;
          if IB.IBP.theta > 0
            IB.Res.Int(Dyzi).rimg = Rimg;
          end
        else
          Pimg = IB.Res.Int(Dyzi).img;
          Yr = IB.Res.Int(Dyzi).Yr;
          Zr = IB.Res.Int(Dyzi).Zr;
          if IB.IBP.theta > 0
            Rimg = IB.Res.Int(Dyzi).rimg;
          end
        end
        if nargout > 0
          figs(Dyzi) = figure;
          h = image(Yr*dyz,Zr*dyz,Pimg','CDataMapping','scaled'); shg
          pos = get(figs(Dyzi),'position');
          pos(1) = 5 + (Dyzi-1)*pos(3);
          pos(2) = 5;
          set(figs(Dyzi),'position',pos);
          set(gca,'DataAspectRatio',[1 1 1],'Ydir','normal');
          ch = colorbar;
          mname = strrep(IB.IBP.mnc,'_','\_');

          title(sprintf('%s: %.1fmm Detector %d passes %d samples', ...
            mname, Dyz*10, IB.IBP.ICOS_passes, ...
            IB.IBP.beam_samples ));
          ylabel(ch,'Nomalized to single ICOS injection');
          xlabel('Y position of detector center, cm');
          ylabel('Z position of detector center, cm');
          drawnow; shg;

          if IB.IBP.theta > 0
            figs(Dyzi+length(IB.IBP.Dyz)) = figure;
            h = image(Yr*dyz,Zr*dyz,Rimg','CDataMapping','scaled'); shg
            pos = get(figs(Dyzi+length(IB.IBP.Dyz)),'position');
            pos(1) = 5 + (Dyzi-1)*pos(3);
            pos(2) = 5 + pos(4);
            set(figs(Dyzi+length(IB.IBP.Dyz)),'position',pos);
            set(gca,'DataAspectRatio',[1 1 1],'Ydir','normal');
            ch = colorbar;
            mname = strrep(IB.IBP.mnc,'_','\_');

            title(sprintf('%s: %.1fmm Detector bad angle', ...
              mname, Dyz*10));
            ylabel(ch,'Nomalized to single ICOS injection');
            xlabel('Y position of detector center, cm');
            ylabel('Z position of detector center, cm');
            drawnow; shg;
          end
        end
      end
      if nargout > 0
        ff = figs;
      end
    end
    
    function analyze_angle(IB)
      % IB.analyze_angle()
      % Reviews the incident angles at the detector and generates a
      % plot showing how the detector power would have been different if
      % the acceptance angle were different. This is somewhat easier to do
      % than running multiple analyses with different target angles.
      
      % I'm taking a shortcut, assuming normal the detector is normal to
      % the optical axis.
      Itheta = acosd(IB.Res.D(1:IB.Res.N,1));
      P = IB.Res.P(1:IB.Res.N);
      Ptotal = sum(P);
      r = max(abs(IB.Res.E(1:IB.Res.N,[2 3])),[],2);
      [Ithetas,I] = sort(Itheta);
      Ps = P(I);
      rs = r(I);
      max_rs = max(rs);
      min_rs = min(rs);
      N = 5;
      rss = linspace(min_rs,max_rs,N);
      figure;
      for i = 1:N
        Inside = rs <= rss(i);
        plot(Ithetas(Inside),100*cumsum(Ps(Inside))/Ptotal);
        hold on;
      end
      hold off;
      legend(num2str(rss','%.3f'),'location','northwest');
      xlabel('Acceptance Angle, degrees');
      ylabel('percent power on detector');
    end
    
    function display(IB)
      % IB.display
      % Overrides how ICOS_beam object is displayed. This is probably not
      % the right way to do this, but it seems to work.
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
    
    function PM_o = draw(IB, varargin)
      % PM = IB.draw(options)
      % Renders the underlying optical model. Options can adjust any of the
      % parameters in the IB.P model definition structure.
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
      if nargout > 0
        PM_o = PM;
      end
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
    
    function savefile(IB)
      % IB.savefile
      % Writes the IB object to a .mat file with the name derived
      % from the mnemonic in IB.IBP.mnc. The filename is prefixed
      % with 'IS_'.
      % See also: ICOS_beam
      fname = sprintf('IB_%s.mat', IB.IBP.mnc);
      save(fname, 'IB');
      fprintf(1, 'ICOS_beam saved to %s\n', fname);
    end
  end
end
