classdef ICOS_search < handle
  % ICOS_search object incorporates the parameter solvers
  % of exparam, autofocus, and ICOS_beam analysis
  properties
    ISP % Starting search parameters
    ISopt % Options
    res1 % Initial ICOS/RIM solutions
    res2 % Autofocus solutions
  end
  methods
    function IS = ICOS_search(varargin)
      IS.ISP.RR1 = []; % If not set, choices are search in ed_rim_catalog
      IS.ISP.Rw1 = [];
      IS.ISP.R1 = [];
      IS.ISP.r1 = [];
      IS.ISP.L = [];
      IS.ISP.R2 = []; % If not set, choices are searched in ispcatalog
      IS.ISopt.mnc = '';
      IS.ISopt.R2_lim = [-inf inf];
      IS.ISopt.RR1_lim = [-inf inf];
      IS.ISopt.L_lim = [0 inf];
      IS.ISopt.RL_lim = [2 inf];
      IS.ISopt.RD1_margin = 1; % cm. When RR1 is set, this is added to Rr1 to determine RD1
      IS.ISopt.D2_margin = 1; % cm. When R2 is set, this is added to r2 to determine D2
      for i=1:2:length(varargin)-1
        fld = varargin{i};
        if isfield(IS.ISP, fld)
          IS.ISP.(fld) = varargin{i+1};
        elseif isfield(IS.ISopt,fld)
          IS.ISopt.(fld) = varargin{i+1};
        else
          error('MATLAB:HUARP:badopt', 'Invalid option: "%s"', fld);
        end
      end
      IS.res1 = [];
      IS.res2 = [];
    end
    
    function search_ICOS_RIM(IS, varargin)
      for i=1:2:length(varargin)-1
        fld = varargin{i};
        if isfield(IS.ISP, fld)
          IS.ISP.(fld) = varargin{i+1};
        else
          error('MATLAB:HUARP:badopt', 'Invalid option: "%s"', fld);
        end
      end
      if isempty(IS.ISP.R2)
        IM2 = ispcatalog;
        IR2 = unique([IM2.R_cm]);
        IR2 = IR2(IR2 > IS.ISopt.R2_lim(1) & IR2 < IS.ISopt.R2_lim(2));
      else
        IM2 = [];
        IR2 = IS.ISP.R2;
      end
      nIR2 = length(IR2);
      
      if isempty(IS.ISP.RR1)
        HM1 = ed_rim_catalog;
        HR1 = unique([HM1.R_cm]);
        HR1 = HR1(HR1 > IS.ISopt.RR1_lim(1) & HR1 < IS.ISopt.RR1_lim(2));
      else
        HM1 = [];
        HR1 = IS.ISP.RR1;
      end
      nHR1 = length(HR1);
      ntrials = nIR2*nHR1;
      
      P.R1 = IS.ISP.R1;
      if ~isempty(IS.ISP.L)
        P.L = IS.ISP.L;
      end
      if ~isempty(IS.ISP.r1)
        P.r1 = IS.ISP.r1;
      end
      P.Rw1 = IS.ISP.Rw1;
      Res = cell(ntrials,0);
      trialn = 1;
      for IR2i = 1:nIR2
        P.R2 = IR2(IR2i);
        for HR1i = 1:nHR1
          P.RR1 = HR1(HR1i);
          Res{trialn} = exparam(P)';
          trialn = trialn+1;
        end
      end
      res = [Res{:}]';
      nres = length(res);
      issane = ones(nres,1);
      if ~isempty(HM1)
        HM1Rcm = [HM1.R_cm];
        HM1Din = [HM1.dia_in];
      end
      if ~isempty(IM2)
        M2Rcm = [IM2.R_cm];
        M2Din = [IM2.dia_in];
      end
      for i = 1:nres
        % Now let's check for some sanity:
        % Is Rr1+0.3 < RD1? (spot radius less than mirror radius)
        % is RR1 < RD1? (radius of curvature less than mirror radius)
        Rr1 = res(i).Rr1;
        RR1 = res(i).RR1;
        if ~isempty(HM1)
          RD1 = HM1Din(HM1Rcm == RR1);
          RD1a = min(RD1(Rr1+0.3 < RD1*2.54/2));
          res(i).RD1 = RD1a;
          if isempty(RD1a)
            fprintf(1,'Discarding %d: Rr1 (%.1f) exceeds RD1/2 (%.1f)\n', i, Rr1, max(RD1)*2.54/2);
            issane(i) = 0;
          end
        else
          res(i).RD1 = (Rr1+IS.ISopt.RD1_margin)*2/2.54;
        end
        
        r2 = res(i).r2;
        R2 = res(i).R2;
        if ~isempty(IM2)
          D2 = M2Din(M2Rcm == R2);
          D2a = min(D2(r2+0.5 < D2*2.54/2));
          res(i).D2 = D2a;
          if isempty(D2)
            fprintf(1,'Discarding %d: r2 (%.1f) exceeds D2/2 (%.1f)\n', i, r2, max(D2)*2.54/2);
            issane(i) = 0;
          end
        else
          res(i).D2 = (r2+IS.ISopt.D2_margin)*2/2.54;
        end
        
        r1 = res(i).r1;
        D1 = 3; % fixed for now
        if r1 > D1*2.54/2
          fprintf(1,'Discarding %d: r1 (%.1f) exceeds D1/2 (%.1f)\n', i, r1, D1*2.54/2);
          issane(i) = 0;
        end
        
        L = res(i).L;
        if L < IS.ISopt.L_lim(1) || L > IS.ISopt.L_lim(2)
          fprintf(1,'Discarding %d: cell length %f out of range\n', i, L);
          issane(i) = 0;
        end
        
        RL = res(i).RL;
        if RL < IS.ISopt.RL_lim(1) || RL > IS.ISopt.RL_lim(2)
          fprintf(1,'Discarding %d: RIM length %f out of range\n', i, RL);
          issane(i) = 0;
        end
      end
      res = res(issane > 0);
      % Need to propagate CT1 and CT2 through the catalog
      i = 0;
      opt_OK = ones(length(res),1)>0;
      while i < length(res)
        %
        i = i+1;
        P = render_model(res(i));
        % fine tune dy/dz
        P.stop_ICOS = 0;
        P.visible = 0;
        P.visibility = 0;
        P.focus = 0;
        P.HR = 0;
        P.ICOS_passes_per_injection = 100;
        P.max_rays = 3000;
        P.plot_endpoints = 0;
        P.evaluate_endpoints = 3;
        P.skip.overlap = 1;
        P.skip.total_power = 1;
        P.skip.mean_angle = 1;
        P.skip.RIM_passes = 1;
        delta = 2.5e-3;
        best_eccentricity = 2;
        eccentricity = 1;
        iteration = 1;
        while delta > 1e-4 && eccentricity < best_eccentricity && iteration < 10
          best_eccentricity = eccentricity;
          PM = ICOS_Model6(P,'dy',P.dy+linspace(-delta,delta,5),'dz',P.dz+linspace(-delta,delta,5));
          PM = PM.clean_results;
          [new_dy,new_dz,eccentricity,ri,ci] = PM.identify_minimum('eccentricity');
          if isempty(new_dy) || isempty(new_dz)
            delta = delta*10;
          elseif eccentricity < best_eccentricity
            P.dy = new_dy;
            P.dz = new_dz;
            if ri ~= 1 && ri ~= size(PM.Results.eccentricity,1) && ...
                ci ~= 1 && ci ~= size(PM.Results.eccentricity,2)
              delta = delta/2;
            end
          end
          iteration = iteration + 1;
        end
        % Now let's see if we can't get the RIM working
        P.stop_ICOS = 1;
        P.HR = 1;
        P.visible = 0;
        P.plot_endpoints = 0;
        P.evaluate_endpoints = 1;
        P.skip.RIM_passes = 0;
        delta = 2; % cm +/-
        iteration = 0;
        while delta > 0.1
          iteration = iteration + 1;
          PM = ICOS_Model6(P,'herriott_spacing', ...
            P.herriott_spacing + linspace(-delta,delta,21));
          criteria = 'eccentricity';
          PM = PM.clean_results;
          PM.Results.eccentricity(PM.Results.RIM_passes <= 1) = NaN;
          [P.herriott_spacing,~,RIM_passes] = PM.identify_minimum('-RIM_passes');
          RIM_passes = abs(RIM_passes)+1;
          [P.herriott_spacing,~,~] = PM.identify_minimum(criteria);
          if isempty(P.herriott_spacing)
            delta = 0;
          else
            delta = delta/10;
          end
        end
        if isempty(P.herriott_spacing)
          fprintf(1,'Solution %d/%d could not optimize herriott_spacing\n', ...
            i, length(res));
          opt_OK(i) = false;
        else
          P.stop_ICOS = 0;
          P.ICOS_passes_per_injection = 20;
          P.max_rays = ceil(RIM_passes+5)*5*P.ICOS_passes_per_injection;
          P.plot_endpoints = 0;
          P.evaluate_endpoints = 3;
          delta = .1;
          PM = ICOS_Model6(P,'herriott_spacing', ...
            P.herriott_spacing + linspace(-delta,delta,11));
          PM = PM.clean_results;
          PM.Results.eccentricity(PM.Results.RIM_passes <= 1) = NaN;
          [P.herriott_spacing,~,~] = PM.identify_minimum(criteria);
          % Calculate overlap
          P.evaluate_endpoints = 2;
          P.skip.total_power = 1;
          P.skip.eccentricity = 1;
          P.skip.mean_angle = 1;
          P.skip.overlap = 0;
          P.skip.RIM_passes = 0;
          P.visible = 0;
          P.HR = 0.98; % restored from before
          PM = ICOS_Model6(P);
          overlap = PM.Results.overlap;
          %
          P.visible = 0;
          P.visibility = [];
          P.focus = 1;
          P.evaluate_endpoints = -1;
          res(i).ORd2 = -P.dy;
          res(i).ORs2 = P.dz;
          res(i).ORL = P.herriott_spacing;
          res(i).overlap = overlap;
          fprintf(1,'Completed %d/%d: RR1:%.1f RL:%.1f r1:%.1f R2:%.1f overlap:%.1f\n', ...
            i, length(res), res(i).RR1, res(i).RL, res(i).r1, res(i).R2, ...
            res(i).overlap);
        end
      end
      IS.res1 = res(opt_OK);
      IS.savefile;
    end
    
    function search_focus(IS, varargin)
      SFopt.select = [];
      for i=1:2:length(varargin)-1
        fld = varargin{i};
        if isfield(IS.ISopt,fld)
          IS.ISopt.(fld) = varargin{i+1};
        elseif isfield(SFopt, fld)
          SFopt.(fld) = varargin{i+1};
        else
          error('MATLAB:HUARP:badopt', 'Invalid option: "%s"', fld);
        end
      end
      if isempty(SFopt.select)
        res = IS.res1;
      else
        res = IS.res1(SFopt.select);
      end
      if ~isempty(IS.res2) && ~isfield(IS.res2,'Nres2')
        for i=1:length(IS.res2)
          IS.res2(i).Nres2 = i;
        end
      end
      i = 0;
      P = ICOS_Model6.props;
      LTS = fields(P.LensTypes);
      results = 0;
      while i < length(res)
        i = i+1;
        P = render_model(res(i), 'visibility', [0 0 0], ...
          'focus', 1, 'ICOS_passes_per_injection', 100, ...
          'max_rays', 3000);
        n = res(i).n; % 2.4361
        d = res(i).d2*n;
        s = res(i).s2;
        r = res(i).r2;
        j = 0;
        Lens1 = '';
        Li = 1;
        if d >= 0
          % place a PM lens at 0.2 that has a large enough diameter
          % and will give us a negative divergence.
          x = 0.2;
          rx = sqrt((r+d.*x).^2 + s.^2.*x.^2);
          dx = (d.*r + (d.^2+s.^2).*x)./rx;
          sx = s.*r./rx;
          f = P.LensTypes.Lens1.EFL;
          dxp = dx - rx/f;
          if dxp < 0
            Lens1 = 'Lens1 @ 0.2 + ';
          else
            Lens1 = 'ERROR ';
          end
          r = rx;
          % d = dxp;
          s = sx;
          P.Lenses = {'Lens1'};
          P.Lens_Space = 0.2;
          Li = 2;
          % recalculate d and continue
          P.HR = 0;
          P.evaluate_endpoints = 5;
          P.visible = 0;
          PM = ICOS_Model6(P);
          [xyz,oxyz] = PM.M.extract_endpoints(5);
          dxyz = xyz-oxyz;
          dxyz = diag(1./dxyz(:,1)) * dxyz;
          ro = sqrt(sum(oxyz(:,2:3).^2,2));
          dyz = sum(dxyz(:,2:3).*oxyz(:,2:3),2)./ro;
          rx = mean(ro);
          d = mean(dyz);
          s = s.*r./rx;
          r = rx;
        end
        while j < length(LTS)
          j = j+1;
          P.Lenses{Li} = LTS{j};
          Lens = P.LensTypes.(LTS{j});
          if strcmp(Lens.type,'positive_meniscus')
            f = Lens.EFL;
          else
            f = -Lens.EFL;
          end
          det_acc_limit = 14.9; % degrees acceptance angle
          [x,theta,ds] = pick_lens_x(r,d,s,f,det_acc_limit,Lens.r);
          if ~isempty(x)
            P.Lens_Space(Li) = x;
            P.detector_spacing = ds;
            P.visible = 0;
            PM = ICOS_Model6(P);
            [xyz,~,div,skew] = PM.M.extract_endpoints_skew(4+Li);
            theta2 = mean(rad2deg(atan(sqrt(div.^2 + skew.^2))));
            if theta2 > det_acc_limit + .01
              ls0 = P.Lens_Space(Li);
              theta0 = theta2;
              delta_ds = 0.4;
              theta1 = theta0;
              while theta1 <= theta0 && theta1 >= det_acc_limit
                ls1 = ls0+delta_ds;
                P.Lens_Space(Li) = ls1;
                PM = ICOS_Model6(P);
                [xyz,~,div,skew] = PM.M.extract_endpoints_skew(4+Li);
                theta1 = mean(rad2deg(atan(sqrt(div.^2 + skew.^2))));
                delta_ds = delta_ds + 0.1;
              end
              theta2 = theta1;
              if theta2 < det_acc_limit
                while abs(theta2-det_acc_limit) > .01
                  ls2 = mean([ls0, ls1]);
                  P.Lens_Space(Li) = ls2;
                  PM = ICOS_Model6(P);
                  [xyz,~,div,skew] = PM.M.extract_endpoints_skew(4+Li);
                  theta2 = mean(rad2deg(atan(sqrt(div.^2 + skew.^2))));
                  if theta2 > det_acc_limit
                    ls0 = ls2;
                  else
                    ls1 = ls2;
                  end
                end
              end
            end
            if theta2 <= det_acc_limit + .01
              div = mean(div);
              skew = mean(skew);
              rx = mean(sqrt(sum(xyz(:,2:3).^2,2)));
              d_ds = -rx*div/(div^2+skew^2);
              ds = ds + d_ds;
              P.detector_spacing = ds;
              if isempty(IS.res2)
                IS.res2 = res(i);
                IS.res2.Nres2 = 1;
                results = 1;
              else
                results = length(IS.res2)+1;
                IS.res2(results).Nres2 = results;
                flds = fieldnames(res);
                for fi=1:length(flds)
                  fld = flds{fi};
                  IS.res2(results).(fld) = res(i).(fld);
                end
              end
              % res_a = res(i);
              IS.res2(results).Lenses = P.Lenses;
              IS.res2(results).Lens_Space = P.Lens_Space;
              IS.res2(results).detector_spacing = ds;
              IS.res2(results).theta = theta2;
              IS.res2(results).r_d = s*r / tand(theta2);
              IS.res2(results).Ltot = IS.res2(results).ORL ...
                + IS.res2(results).L + sum(IS.res2(results).Lens_Space) ...
                + IS.res2(results).detector_spacing;
              fprintf(1,'%d: (%d,%d) %s%s @ %.2f (%.1f deg)\n', results, i, j, Lens1, LTS{j}, x, theta);
            else
              warning('MATLAB:HUARP:BadAngle', 'Rejecting solution for bad angle');
            end
          end
        end
      end
      IS.savefile;
    end
    
    function analyze(IS, varargin)
      % IS.analyze(opts);
      % opts include:
      % 'ICOS_passes', n
      % 'Nsample', n
      % 'opt_n', n
      % plus and ICOS_beam options.
      if isempty(IS.res2)
        error('MATLAB:HUARP:NoFocus', ...
          'Must invoke IS.search_focus() before IS.analyze()');
      end
      Opt.ICOS_passes = 50;
      Opt.Nsamples = 100;
      % Opt.Herriott_passes = 100;
      % Opt.n_optics = 5;
      % Opt.rng_state = [];
      Opt.opt_n = [];
      Opt.select = [];
      IBopt = {};
      for i=1:2:length(varargin)-1
        fld = varargin{i};
        if isfield(IS.ISopt,fld)
          IS.ISopt.(fld) = varargin{i+1};
        elseif isfield(Opt, varargin{i})
          Opt.(varargin{i}) = varargin{i+1};
        else
          IBopt{end+1} = varargin{i};
          IBopt{end+1} = varargin{i+1};
        end
      end
      if ~isempty(IS.res2) && ~isfield(IS.res2,'Nres2')
        for i=1:length(IS.res2)
          IS.res2(i).Nres2 = i;
        end
      end
      % Use unique to map between IS.res2.Nres2 and indices
      [NR2,INR2] = unique([IS.res2.Nres2]);
      if isempty(Opt.select)
        Opt.select = NR2;
      end
      NRi = interp1(NR2,INR2,Opt.select,'nearest');
      NRi = NRi(~isnan(NRi));
      ff = [];
      if iscolumn(NRi)
        NRi = NRi';
      end
      for i=NRi
        ofile = sprintf('IS_%s.%d_%dx%d', IS.ISopt.mnc, ...
          IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
        P = render_model(IS.res2(i));
        n_optics = 4 + length(P.Lenses);
        Track_Power = 1;
        if ~isempty(Opt.opt_n)
          opt_n = Opt.opt_n;
          Track_Power = 0;
          ofile = sprintf('%s_%d', ofile, opt_n);
        else
          opt_n = n_optics;
        end
        ofile = sprintf('%s.mat', ofile);
        if ~exist(ofile, 'file')
          P.visible = 0;
          % P.HR = 0;
          P.focus = 1;
          P.evaluate_endpoints = -1;
          IB = ICOS_beam(@ICOS_Model6, P);
          IB.Sample('beam_samples', Opt.Nsamples, ...
            'ICOS_passes', Opt.ICOS_passes, 'opt_n', opt_n, ...
            'n_optics', n_optics, 'Track_Power', Track_Power, IBopt{:});
          if opt_n == n_optics
            ff2 = IB.Integrate;
            if ~isempty(ff)
              delete(ff);
              drawnow;
            end
            ff = ff2;
          end
          save(ofile, 'IB');
          fprintf(1, 'ICOS_beam %d Saved result to %s\n', i, ofile);
          Pwr = IB.PowerSummary;
          IS.res2(i).NH = Pwr.NH;
          if isfield(Pwr,'max_pwr')
            IS.res2(i).max_pwr = Pwr.max_pwr;
          elseif ~isfield(IS.res2,'max_pwr')
            IS.res2(i).max_pwr = [];
          end
          IS.savefile;
        else
          fprintf(1, 'Skipping result %d\n', i);
        end
      end
    end
    
    function savefile(IS)
      fname = sprintf('IS_%s.mat', IS.ISopt.mnc);
      save(fname, 'IS');
      fprintf(1, 'ICOS_search saved to %s\n', fname);
    end
  end
end