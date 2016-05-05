classdef ICOS_search < handle
  % ICOS_search object incorporates the parameter solvers
  % of exparam, autofocus, and ICOS_beam analysis. The basic concept is to
  % search through standard elements from an optical catalog to find useful
  % configurations.
  %
  % See also:
  %  ICOS_search.ICOS_search
  %  ICOS_search.search_ICOS_RIM
  %  ICOS_search.search_focus (deprecated)
  %  ICOS_search.search_focus2
  %  ICOS_search.explore_focus
  %  ICOS_search.analyze
  %  ICOS_search.savefile
  properties
    ISP % Starting search parameters
    ISopt % Options
    res1 % Initial ICOS/RIM solutions
    res2 % Autofocus solutions
  end
  methods
    function IS = ICOS_search(varargin)
      % IS = ICOS_search(options)
      % The options define the search parameters via keyword/value
      % pairs.
      %   mnc: required, used for filenames
      %   Configuration constraints: RR1, Rw1, R1, r1, L, R2
      %     5 of these need to be known to unambiguously define
      %     an ICOS configuration. RR1 and R2, if not specified,
      %     can be obtained by trying elements from their
      %     respective catalogs. I believe this program currently requires
      %     Rw1 and R1 be defined plus either r1 or L. Current research is
      %     looking at defining other sets.
      %   ICOS_search options: These control the catalog searches:
      %     R2_lim: [min max] defines the range of values for R2
      %     RR1_lim: [min max] same for RR1
      %     L_lim: [min max] range for L
      %     RL_lim: [min max] range for RL
      %     RD1_margin: Additional radial distance beyond Rr1 required when
      %          determining RD1. Defaults to 1 cm.
      %     D2_margin: Additional radial distance beyond r2 required when
      %          determining D2. Defaults to 1 cm.
      %     focus_visible: Used in search_focus2. If 0, all visualizations
      %          are suppressed. If 1, final focus visualizations are
      %          shown. If 2, visualizations are displayed at each step.
      IS.ISP.RR1 = []; % If not set, choices are search in ed_rim_catalog
      IS.ISP.Rw1 = [];
      IS.ISP.R1 = [];
      IS.ISP.r1 = [];
      IS.ISP.L = [];
      IS.ISP.R2 = []; % If not set, choices are searched in ispcatalog
      IS.ISP.optics_n = []; % index of refraction
      IS.ISP.mirrors_n = [];
      IS.ISP.lenses_n = [];
      IS.ISopt.mnc = '';
      IS.ISopt.R2_lim = [-inf inf];
      IS.ISopt.RR1_lim = [-inf inf];
      IS.ISopt.L_lim = [0 inf];
      IS.ISopt.RL_lim = [2 inf];
      IS.ISopt.RD1_margin = 2; % cm. When RR1 is set, this is added to Rr1 to determine RD1
      IS.ISopt.RD1_lip = 0;
      IS.ISopt.RD1_resolution = 2.54/4;
      IS.ISopt.D1_margin = 0.3; % cm. When R1 is set, this is added to r1 to determine D1
      IS.ISopt.D1_lip = 0.125*2.54; % cm on the radius.
      IS.ISopt.D1_resolution = 2.54/4; % Applies to radius
      IS.ISopt.D2_margin = 0.3; % cm. When R2 is set, this is added to r2 to determine D2
      IS.ISopt.D2_lip = 0.125*2.54; % cm on the radius
      IS.ISopt.D2_resolution = 2.54/4; % applies to radius
      IS.ISopt.focus_visible = 1; % 0: don't draw any focus 1: just finished focus 2: each iteration
      IS.ISopt.max_focus_length = 20; % abandon is sum of lens space exceeds
      IS.ISopt.allow_negative_focus = 0;
      IS.ISopt.allow_nondecreasing_focus = 0; % Allow big focus lenses after small ones
      IS.ISopt.beam_diameter = 0.4;
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
      % IS.search_ICOS_RIM(options)
      % Options can include any of the ICOS configuration options
      % from ICOS_search.ICOS_search.
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
      if ~isempty(IS.ISP.mirrors_n)
        P.n = IS.ISP.mirrors_n;
      elseif ~isempty(IS.ISP.optics_n)
        P.n = IS.ISP.optics_n;
      end
      if ~isempty(IS.ISP.Rw1)
        P.Rw1 = IS.ISP.Rw1;
      end
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
          res(i).RD1 = ...
            (ceil((Rr1+IS.ISopt.RD1_margin+IS.ISopt.RD1_lip)/ ...
              IS.ISopt.RD1_resolution)*IS.ISopt.RD1_resolution - ...
              IS.ISopt.RD1_lip) * 2/2.54;
        end
        
        res(i).D1 = (ceil((res(i).r1+IS.ISopt.D1_margin+IS.ISopt.D1_lip)/...
          IS.ISopt.D1_resolution) * IS.ISopt.D1_resolution ...
          - IS.ISopt.D1_lip) * 2 / 2.54;
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
          % res(i).D2 = (r2+IS.ISopt.D2_margin)*2/2.54;
          res(i).D2 = (ceil((res(i).r2+IS.ISopt.D2_margin+IS.ISopt.D2_lip)/...
            IS.ISopt.D2_resolution) * IS.ISopt.D2_resolution ...
            - IS.ISopt.D2_lip) * 2 / 2.54;
        end
        
        r1 = res(i).r1;
        D1 = res(i).D1;
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
        P = render_model(res(i),'beam_diameter',IS.ISopt.beam_diameter);
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
          if isempty(P.herriott_spacing)
            fprintf(1,'Solution %d/%d could not optimize herriott_spacing with ICOS\n', ...
              i, length(res));
            opt_OK(i) = false;
          else
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
      end
      IS.res1 = res(opt_OK);
      IS.savefile;
    end
    
    function search_focus2(IS, varargin)
      % IS.search_focus2(options)
      % The main relevant options are:
      %   select: takes a list of indices into the IS.res1 array to
      %     be analyzed.
      %   focus_visible: 0 suppresses plots, 1 plots accepted focuses,
      %     2 shows all attempts.
      %   max_lenses: defaults to 3
      %   fix_lenses: cell array of specific lenses to use
      %   injection_scale: Used to scale up radius
      %
      % search_focus2 attempts to build configurations using lenses
      % defined in the ICOS_Model6.props LensTypes array.
      function optimize_focus(IS, resn, P, d, s, r, th, dth, depth, fix, max_lens_radius)
        % optimize_focus(IS, resn, P, d, s, r, th, dth, fix, max_lens_radius)
        % IS ICOS_search object
        % resn The res1 index we are working on
        % P ICOS_Model6 properties
        % d divergence at current last optic
        % s skew at current last optic
        % r beam radius at current last optic
        % th target angle
        % dth target angle tolerance
        % depth is a limit on how many lenses we should allow
        % fix cell array if non-empty, defines specific lenses to use
        % max_lens_radius
        % Assumption is that we need to be converging.
        if th < 0 || th > 90
          error('MATLAB:HUARP:badangle', ...
            'Angle th (%f) must be between 0 and 90', th);
        end
        if dth <= 0
          error('MATLAB:HUARP:badangle', ...
            'Angle dth (%f) must be >= 0', dth);
        end
        if sum(P.Lens_Space) + P.detector_spacing > IS.ISopt.max_focus_length
          fprintf(1,'Abandoning focus for excessive length\n');
          return;
        end
        theta = atand(sqrt(d^2+s^2));
        if d < 0 && theta > th-dth && theta < th+dth
          result = length(IS.res2)+1;
          if isempty(IS.res2)
            IS.res2 = IS.res1(resn);
            IS.res2.Nres2 = 1;
          else
            IS.res2(result).Nres2 = result;
            flds = fieldnames(IS.res1);
            for fi=1:length(flds)
              fldnm = flds{fi};
              IS.res2(result).(fldnm) = IS.res1(resn).(fldnm);
            end
          end
          % res_a = res(i);
          IS.res2(result).Lenses = P.Lenses;
          IS.res2(result).Lens_Space = P.Lens_Space;
          IS.res2(result).detector_spacing = P.detector_spacing;
          IS.res2(result).theta = theta;
          IS.res2(result).r_d = s*r / tand(theta);
          IS.res2(result).Ltot = IS.res2(result).ORL ...
            + IS.res2(result).L + sum(IS.res2(result).Lens_Space) ...
            + IS.res2(result).detector_spacing;
          IS.res2(result).sel = 0;
          fprintf(1,'%d: (%d) Focused: theta=%.2f', result, resn, theta);
          for Li=1:length(P.Lenses)
            fprintf(1, ' %s', P.Lenses{Li});
          end
          fprintf(1, '\n');
          if IS.ISopt.focus_visible > 0
            P.visible = 1;
          else
            P.visible = 0;
          end
          PM = ICOS_Model6(P);
          if P.visible
            title(sprintf('Focus %d', result));
            drawnow;
          end
          return;
        end
        if depth <= 0
          fprintf(1,'Abandoning focus after %d lenses\n', ...
            length(P.Lenses));
          return; % No results
        end
        n_lenses = length(P.Lenses)+1;
        if ~isempty(fix)
          if n_lenses <= length(fix)
            LTS = { fix{n_lenses} };
          else
            LTS = {};
          end
        else
          LTS = fields(P.LensTypes);
        end
        if n_lenses == 1
          if max_lens_radius > 0
            Prev_lens_r = max_lens_radius;
          else
            Prev_lens_r = P.r2;
          end
        else
          Prev_lens_r = P.LensTypes.(P.Lenses{n_lenses-1}).r;
        end
        for LTSi = 1:length(LTS)
          P.Lenses{n_lenses} = LTS{LTSi};
          Lens = P.LensTypes.(LTS{LTSi});
          f = Lens.EFL;
          if strcmp(Lens.type,'negative_meniscus') && f > 0
            f = -f;
            warning('MATLAB:HUARP:negnotneg', ...
              'Negative meniscus %s has positive EFL', LTS{LTSi});
          end
          if (IS.ISopt.allow_negative_focus || f > 0) && ...
              (max_lens_radius == 0 || ...
                Lens.r <= max_lens_radius) && ...
              (IS.ISopt.allow_nondecreasing_focus || ...
                Lens.r <= Prev_lens_r)
            [x,~,ds] = pick_lens_x2(r,d,s,f,th,Lens.r-0.3);
            % pick_lens_x2 could have found
            %   a: no solution -- just try the next lens
            %   b: an incomplete solution (at xmin or xmax)
            %   c: a good solution
            % in cases b and c, we need to verify the angle in the
            % full model. In b, we do not expect a good value, but
            % with c, it is possible we might be a little bit off.
            % In that case, we should try to optimize by binary
            % search.
            if ~isempty(x)
              optimize_rays(IS, P, x, d, ds, th, dth, f, resn, depth, fix, max_lens_radius);
              if length(x) == 3 && x(1) > 5
                % Try for a shorter focus by using minimum
                optimize_rays(IS, P, x(2), d, ds, th, dth, f, resn, depth, fix, max_lens_radius);
              end
            end
          end
        end
      end
      
      function optimize_rays(IS, P, x, d, ds, th, dth, f, resn, depth, fix, max_lens_radius)
        n_lenses = length(P.Lenses);
        P.Lens_Space(n_lenses) = x(1);
        P.detector_spacing = ds;
        if IS.ISopt.focus_visible > 1
          P.visible = 1;
        else
          P.visible = 0;
        end
        P.view = [0 0];
        P.evaluate_endpoints = -1;
        PM = ICOS_Model6(P);
        [oxyz,r2,div,skew] = PM.M.extract_origin_skew(4+n_lenses);
        if isempty(div)
          return;
        end
        theta2 = atand(sqrt(mean(div).^2 + mean(skew).^2));
        % Now do we need to optimize further? Only if
        % length(x) == 3 and theta2 is not close to th
        if length(x) == 3 && (theta2 < th-dth || theta2 > th+dth)
          xmin = x(2);
          xmax = x(3);
          xminok = false;
          xmaxok = false;
          thetamax = -1;
          thetamin = -1;
          while theta2 < th-dth || theta2 > th+dth
            fprintf(1,'Optimizing: f=%f theta=%.2f x=[%.2f %.2f %.2f]\n', ...
              f, theta2, xmin, P.Lens_Space(n_lenses), xmax);
            if theta2 < 0 % try to backtrack toward an 'OK' limit
              if xminok
                xmax = P.Lens_Space(n_lenses);
                xmaxok = false;
                fprintf(1,'Retreating toward xmin\n');
              elseif xmaxok
                xmin = P.Lens_Space(n_lenses);
                xminok = false;
                fprintf(1,'Retreating toward xmax\n');
              else
                fprintf(1,'Abandoning this line: theta<0 and no ok limits\n');
                return;
              end
              % elseif (f > 0 && theta2 < th-dth) || (f < 0 && theta2 > th+dth)
            elseif xor(d > 0, xor(f > 0, theta2 > th+dth))
              % move in
              xmax = P.Lens_Space(n_lenses);
              xmaxok = true;
              thetamax = theta2;
            else
              % move out
              xmin = P.Lens_Space(n_lenses);
              xminok = true;
              thetamin = theta2;
            end
            if xmaxok && xminok
              xnew = interp1([thetamin thetamax],[xmin xmax],th);
            else
              xnew = mean([xmax, xmin]);
            end
            if abs(xnew-P.Lens_Space(n_lenses)) < .01
              fprintf(1,'Abandoning this line: dx < .01\n');
              return;
            end
            P.Lens_Space(n_lenses) = xnew;
            PM = ICOS_Model6(P);
            [oxyz,r2,div,skew] = PM.M.extract_origin_skew(4+n_lenses);
            div = mean(div);
            if div >= 0
              theta2 = -1;
            else
              theta2 = atand(sqrt(div.^2 + mean(skew).^2));
            end
          end
        end
        % [oxyz,r2,div,skew] = PM.M.extract_origin_skew(4+n_lenses);
        r2 = mean(r2);
        d2 = mean(div);
        s2 = mean(skew);
        P.detector_spacing = mean(oxyz(:,1)) - PM.M.Optic{3+n_lenses}.O(1) - ...
          PM.M.Optic{3+n_lenses}.CT - r2*d2/(d2^2+s2^2);
        optimize_focus(IS, resn, P, d2, s2, r2, th, dth, depth-1, fix, max_lens_radius);
      end
      
      SFopt.select = [];
      SFopt.det_acc_limit = 14.9;
      SFopt.det_acc_limit_tolerance = 0.05;
      SFopt.max_lenses = 3;
      SFopt.fix_lenses = {};
      SFopt.injection_scale = 1;
      SFopt.max_lens_radius = 0;
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
      % P = ICOS_Model6.props;
      results = 0;
      while i < length(res)
        i = i+1;
        P = render_model(res(i), 'visibility', [0 0 0], ...
          'focus', 1, 'ICOS_passes_per_injection', 100, ...
          'injection_scale', SFopt.injection_scale, ...
          'max_rays', 3000, 'HR', 0, ...
          'beam_diameter',IS.ISopt.beam_diameter);
        if ~isempty(IS.ISP.lenses_n)
          P.lenses_n = IS.ISP.lenses_n;
        elseif ~isempty(IS.ISP.optics_n)
          P.lenses_n = IS.ISP.optics_n;
        end
        n = res(i).n; % 2.4361
        d = res(i).d2*n;
        s = res(i).s2;
        r = res(i).r2;
        optimize_focus(IS, i, P, d, s, r, SFopt.det_acc_limit, ...
          SFopt.det_acc_limit_tolerance, SFopt.max_lenses, ...
          SFopt.fix_lenses, SFopt.max_lens_radius);
      end
      IS.savefile;
    end
    
    function search_focus(IS, varargin)
      % IS.search_focus(options) (deprecated)
      % The main relevant option is 'select', which takes a list of
      % indices into the IS.res1 array to be analyzed.
      %
      % search_focus attempts to build configurations using lenses
      % defined in the ICOS_Model6.props LensTypes array. It does
      % not optimize the result very well, which is why search_focus2 was
      % created.
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
      % P = ICOS_Model6.props;
      LTS = fields(P.LensTypes);
      results = 0;
      while i < length(res)
        i = i+1;
        P = render_model(res(i), 'visibility', [0 0 0], ...
          'focus', 1, 'ICOS_passes_per_injection', 100, ...
          'injection_scale', SFopt.injection_scale, ...
          'max_rays', 3000, ...
          'beam_diameter',IS.ISopt.beam_diameter);
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
    
    function explore_focus(IS, plotnum)
      if isempty(IS.res2)
        error('MATLAB:HUARP:NoFocus', ...
          'Must invoke IS.search_focus2() before IS.explore_focus()');
      end
      if ~isfield(IS.res2,'sel')
        for i=1:length(IS.res2)
          IS.res2(i).sel = 0;
        end
      end
      DS = [IS.res2.detector_spacing];
      LS = zeros(size(DS));
      for i=1:length(LS)
        LS(i) = IS.res2(i).Lens_Space(end);
      end
      Lpos = LS ./ (LS+DS);
      Ltot = [IS.res2.Ltot];
      
      L1EFL = ones(length(IS.res2),1);
      L2r = ones(length(IS.res2),1);
      P = ICOS_Model6.props;
      for i=1:length(IS.res2)
        Ltype = IS.res2(i).Lenses{1};
        L1EFL(i) = P.LensTypes.(Ltype).EFL;
        Ltype = IS.res2(i).Lenses{end};
        L2r(i) = P.LensTypes.(Ltype).r;
      end
      
      sel = [IS.res2.sel];
      v = sel ~= 0;
      
      f = figure;
      if nargin < 2 || plotnum == 1
        h = plot(Ltot,Lpos,'.',Ltot(v),Lpos(v),'or');
        xlabel('L_{total} cm');
        ylabel('Relative lens position');
        title(sprintf('%s Focus Analysis', strrep(IS.ISopt.mnc,'_','\_')));
        hdt = datacursormode;
        set(hdt,'UpdateFcn', ...
          {@ICOS_search.data_cursor_text_func,Ltot,Lpos,Ltot,Lpos,L1EFL,L2r});
        set(h,'buttondownfcn', @(s,e) ICOS_search.ex_bdf(s,e,IS,Ltot,Lpos));
        set(gca, 'buttondownfcn', @(s,e) ICOS_search.ex_bdf(s,e,IS,Ltot,Lpos));
      elseif plotnum == 2
      end
      waitfor(f);
    end
    
    function analyze(IS, varargin)
      % IS.analyze(options);
      % options include:
      %   'ICOS_passes', n
      %   'Nsamples', n
      %   'HR', Herriott mirror reflectivity (0.98)
      %   'T', ICOS mirror transmission (250ppm)
      %   'opt_n', n
      %   'select', array of Nres2 indexes to analyze
      % plus any ICOS_beam options
      %
      % IS.analyze runs the ICOS_beam Sample and Integrate analyses for all
      % the specified configurations and stores the result within the
      % IS.res2 array. The IS object is saved after each ICOS_beam
      % analysis, so the process can be interrupted without losing too much
      % work. On restart, the existence of saved ICOS_beam objects will
      % allow this process to skip ahead.
      if isempty(IS.res2)
        error('MATLAB:HUARP:NoFocus', ...
          'Must invoke IS.search_focus() before IS.analyze()');
      end
      Opt.ICOS_passes = 50;
      Opt.Nsamples = 100;
      Opt.HR = 0.98;
      Opt.T = 250e-6;
      Opt.Herriott_passes = 100;
      Opt.opt_n = [];
      Opt.select = [];
      if isfield(IS.res2,'sel')
        seli = find([IS.res2.sel]);
        for i=1:length(seli)
          seli(i) = IS.res2(seli(i)).Nres2;
        end
        Opt.select = seli;
      end
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
      if length(NR2) > 1
        NRi = interp1(NR2,INR2,Opt.select,'nearest');
        NRi = NRi(~isnan(NRi));
      elseif ~isempty(NR2) && Opt.select(1) == NR2
        Opt.select = NR2;
        NRi = INR2;
      else
        Opt.select = [];
        NRi = [];
      end
      ff = [];
      if iscolumn(NRi)
        NRi = NRi';
      end
      for i=NRi
        IBmnc = sprintf('%s.%d_%dx%d', IS.ISopt.mnc, ...
          IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
%         ofile = sprintf('IB_%s.%d_%dx%d', IS.ISopt.mnc, ...
%           IS.res2(i).Nres2, Opt.ICOS_passes,Opt.Nsamples);
        P = render_model(IS.res2(i),'beam_diameter',IS.ISopt.beam_diameter);
        if ~isempty(IS.ISP.lenses_n)
          P.lenses_n = IS.ISP.lenses_n;
        elseif ~isempty(IS.ISP.optics_n)
          P.lenses_n = IS.ISP.optics_n;
        end
        P.HR = Opt.HR;
        if Opt.HR == 0
          Opt.Herriott_passes = 1;
        end
        P.T = Opt.T;
        n_optics = 4 + length(P.Lenses);
        Track_Power = 1;
        if ~isempty(Opt.opt_n)
          opt_n = Opt.opt_n;
          Track_Power = 0;
          IBmnc = sprintf('%s_%d', IBmnc, opt_n);
        else
          opt_n = n_optics;
        end
        ofile = sprintf('IB_%s.mat', IBmnc);
        if ~exist(ofile, 'file')
          P.visible = 0;
          % P.HR = 0;
          P.focus = 1;
          P.evaluate_endpoints = -1;
          IB = ICOS_beam(@ICOS_Model6, P);
          IB.Sample('beam_samples', Opt.Nsamples, ...
            'ICOS_passes', Opt.ICOS_passes, 'opt_n', opt_n, ...
            'n_optics', n_optics, 'Track_Power', Track_Power, ...
            'Herriott_passes', Opt.Herriott_passes, ...
            'mnc', IBmnc, IBopt{:});
          if opt_n == n_optics
            IB.Integrate;
%             ff2 = IB.Integrate;
%             if ~isempty(ff)
%               delete(ff);
%               drawnow;
%             end
%             ff = ff2;
          end
          IB.savefile;
          % save(ofile, 'IB');
          % fprintf(1, 'ICOS_beam %d Saved result to %s\n', i, ofile);
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
      % IS.savefile
      % Writes the IS object to a .mat file with the name derived
      % from the mnemonic in IS.ISopt.mnc. The filename is prefixed
      % with 'IS_'.
      % See also: ICOS_beam
      fname = sprintf('IS_%s.mat', IS.ISopt.mnc);
      save(fname, 'IS');
      fprintf(1, 'ICOS_search saved to %s\n', fname);
    end
  end
  
  methods(Static)
    function output_txt = data_cursor_text_func(~,event_obj,...
        x,y,Ltot,Lpos,L1EFL,L2r)
      % Display the position of the data cursor
      % obj          Currently not used (empty)
      % event_obj    Handle to event object
      % output_txt   Data cursor text string (string or cell array of strings).
      
      pos = get(event_obj,'Position');
      % output_txt = {['X: ',num2str(pos(1),4)],...
      %    ['Y: ',num2str(pos(2),4)]};

      % If there is a Z-coordinate in the position, display it as well
      % if length(pos) > 2
      %   output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
      % end
      
      output_txt = {};
      i = find(pos(1)==x & pos(2)==y,1);
      if isempty(i)
        output_txt = {'<no text>'};
      else
        output_txt{end+1} = sprintf('Ltot = %.2f cm', Ltot(i));
        output_txt{end+1} = sprintf('Lpos = %.2f', Lpos(i));
        output_txt{end+1} = sprintf('L1EFL = %.2f cm', L1EFL(i));
        output_txt{end+1} = sprintf('L2r = %.2f cm', L2r(i));
      end
    end
    
    function ex_bdf(src,~,IS,X,Y)
      % Button Down Function for explore scripts
      srctype = get(src,'type');
      while ~strcmp(srctype,'axes')
        src = get(src,'parent');
        if isempty(src)
          fprintf(1,'No axes found\n');
          return;
        end
        srctype = get(src,'type');
      end
      set(gcf,'units','normalized');
      axesHandle  = get(src,'Parent');
      coordinates = get(axesHandle,'CurrentPoint');
      coordinates = coordinates(1,1:2);
      % fprintf(1,'Click at coordinates: %f, %f\n', coordinates(1), coordinates(2));
      rect = rbbox([coordinates 0 0]);
      % fprintf(1,'Final rect: [%.2f %.2f %.2f %.2f]\n', rect);
      xl = get(gca,'xlim');
      yl = get(gca,'ylim');
      ap = get(gca,'position');
      XL = interp1([ap(1) ap(1)+ap(3)],xl,[rect(1),rect(1)+rect(3)],'linear','extrap');
      YL = interp1([ap(2) ap(2)+ap(4)],yl,[rect(2),rect(2)+rect(4)],'linear','extrap');
      % fprintf(1,'Final rect: X: [%.2f, %.2f] Y: [%.2f, %.2f]\n', ...
      %   XL, YL);
      v = find(X >= XL(1) & X <= XL(2) & Y >= YL(1) & Y <= YL(2));
      for i=1:length(v)
        IS.res2(v(i)).sel = ~IS.res2(v(i)).sel;
      end
      sel = find([IS.res2.sel]);
      h = get(gca,'children');
      if length(h) == 2
        m = get(h,'marker');
        hi = find(strcmp(m, 'o'));
        set(h(hi),'XData',X(sel),'YData',Y(sel));
      else
        hold on;
        plot(X(sel),Y(sel),'or');
        hold off;
      end
    end
  end
end
