classdef ICOS_Model6 < opt_model_p
  % This model will include a Herriot Mirror
  % Optic:
  %  1: Herriott Mirror
  %  2: HR mirror 1
  %  3: HR mirror 2
  %  *: Lenses!
  %  x: Detector
  % ICOS_Model6 supports different radii of curvature for the
  % two ICOS mirrors.
  properties
  end
  
  methods
    function PM = ICOS_Model6(P,varargin)
      PM = PM@opt_model_p(P,varargin{:});
    end
    
    function drawn = visualize(PM, P)
      drawn = false;
      if PM.M.visible
        drawn = PM.visualize@opt_model_p(P);
      elseif P.plot_endpoints > 0
        PM.M.plot_endpoints(P.plot_endpoints);
        PM.draw_iteration_title(P);
        drawn = true;
      end
    end
    
    function Res = evaluate_endpoints(PM,P)
      Res = PM.results_struct;
      if P.evaluate_endpoints < 0
        return;
      elseif P.evaluate_endpoints > 0
        opt_n = P.evaluate_endpoints;
      elseif P.plot_endpoints > 0
        opt_n = P.plot_endpoints;
      else
        opt_n = length(PM.M.Optic);
      end

      if ~(P.skip.overlap && P.skip.eccentricity)
        xyz = PM.M.extract_endpoints(opt_n);
        Res.n_rays = size(xyz,1);
      end
      
      Res = PM.results_struct;
      Res.inside = PM.M.Inside;

      % Total power calculation
      if ~P.skip.total_power
        vf = find([PM.M.Rays(1:PM.M.n_rays).n_opt] == opt_n);
        Res.total_power = 0;
        for i = 1:length(vf)
          ray = PM.M.Rays(vf(i)).ray;
          if ray.Inside
            Res.total_power = Res.total_power + ray.P;
          end
        end
      end

      % Overlap calculation
      % Now ignores opt_n and always analyzes for optic #2,
      % the first ICOS mirror.
      if ~P.skip.overlap
        n_opt = [PM.M.Rays(1:PM.M.n_rays).n_opt];
        v = n_opt == 2; % was opt_n
        vi = find(v);
        RIM_pass = zeros(length(vi),1);
        for i=1:length(vi)
          RIM_pass(i) = PM.M.Rays(vi(i)).ray.pass(1);
        end
        RIM_passes = unique(RIM_pass);
        overlap = zeros(size(RIM_passes));
        for j=1:length(RIM_passes)
          ji = find(RIM_pass == RIM_passes(j));
          n_pts = length(ji);
          if n_pts > P.n_overlap_spots
            n_pts = P.n_overlap_spots;
            ji = ji(1:n_pts);
          end
          if n_pts > 1
            P0 = PM.M.Rays(vi(ji(1))).ray.E(1,2:3);
            for i=2:n_pts
              d = PM.M.Rays(vi(ji(i))).ray.E(1,2:3) - P0;
              d = sqrt(sum(d.^2,2));
              overlap(j) = overlap(j) + ...
                sum(max(0,P.beam_diameter-d))/P.beam_diameter;
            end
          end
        end
        Res.overlap = mean(overlap);
      end
      
      if ~P.skip.eccentricity
        yz = xyz(:,[2 3]);
        rr = minmax(sqrt(sum(yz.^2,2))');
        Res.max_radius = rr(2);
        Res.eccentricity = sqrt(1 - (rr(1)^2 / rr(2)^2));
      end
      
      if ~P.skip.mean_angle
        Res.mean_angle = mean(PM.M.extract_angles(opt_n));
        if Res.mean_angle > pi
          Res.mean_angle = Res.mean_angle - 2*pi;
        end
        Res.mean_angle = rad2deg(Res.mean_angle);
      end
      
      % RIM_passes
      if ~P.skip.RIM_passes
        if P.HR > 0
          n_inc = [PM.M.Rays(2:PM.M.n_rays).n_inc];
          s_opt = [PM.M.Rays(n_inc).n_opt]; % source optic
          Res.RIM_passes = sum(s_opt == 1);
        else
          Res.RIM_passes = 0;
        end
      end
      
      Res.max_rays = PM.M.max_rays;
    end
    
    % This is a candidate for promotion to opt_model_p:
    function [x, y, z, ri, ci] = identify_minimum(PM,metric)
      if metric(1) == '-'
        invert = -1;
        metric = metric(2:end);
      else
        invert = 1;
      end
      if ~isfield(PM.Results,metric)
        error('MATLAB:HUARP:NotAMetric', ...
          '"%s" is not a metric in this model', metric);
      end
      if isempty(PM.Results.x)
        error('MATLAB:HUARP:NoIndepententVariables', ...
          'There are no independent variables in this model');
      end
      metvals = invert * PM.Results.(metric);
      z = nanmin(nanmin(metvals));
      [i,j] = find(metvals == z,1);
      x = PM.Results.(PM.Results.x)(i,j);
      if ~isempty(PM.Results.y)
        y = PM.Results.(PM.Results.y)(i,j);
      else
        y = [];
      end
      if nargout > 3
        ri = j;
      end
      if nargout > 4
        ci = i;
      end
    end
  end
  
  methods (Static)
    function P = props
      P.optics_n = 2.4361;
      % Herriott Mirror
      P.herriott_spacing = 10; % Before the first ICOS mirror
      P.HRC = 15*2.54; % Herriott radius of curvature
      P.HCT = 0.2;
      P.Hr = 1.5*2.54; % Herriott radius
      P.HR = 0.98; % Herriott reflectivity
      P.C = 3000; % Laser coherence length or equivalent, cm
      
      % ICOS Mirrors
      P.T = 250e-6; % ICOS mirror transmission
      P.mirror_spacing = 50; % mirror spacing
      P.r1 = 1.5*2.54; % mirror radius
      P.R1 = 75;
      P.CT1 = 0.2; % mirror center thickness
      P.r2 = .75*2.54;
      P.R2 = -25.4;
      P.CT2 = .22;

      % Focusing lenses
      P.Lenses = {};
      P.Lens_Space = [];
      
      % Detector
      P.detector_spacing = 3;
      P.D_l = .2; % detector edge size
      P.D_dY = 0; % detector horizontal displacement
      
      % Injection angle
      P.injection_scale = 1; % Scales y0, dy and dz
      % P.y0,z0 are the y and z position of the
      % entrance beam at the back of the first
      % ICOS mirror.
      P.y0 = P.Hr-0.54; % Location of Herriott hole
      P.z0 = 0;
      % P.beam_dy,dz are positional offsets from y0,z0
      % used by ICOS_beam to simulate a beam of finite
      % width;
      P.beam_dy = 0;
      P.beam_dz = 0;
      % P.dy,dz are the direction components of the
      % entrance beam that intersects the first
      % ICOS mirror at P.y0,z0. Assumint P.z0=0,
      % P.dy is the divergence and P.dz is the skew.
      P.dy = .04;
      P.dz = .03;
      
      % Lens types:
      % Custom 3" poitive meniscus
      P.LensTypes.Lens1.type = 'positive_meniscus';
      P.LensTypes.Lens1.r = 3*2.54/2;
      P.LensTypes.Lens1.R1 = 8.0122;
      P.LensTypes.Lens1.R2 = 29.8275;
      P.LensTypes.Lens1.CT = 0.9;
      P.LensTypes.Lens1.EFL = 7.62;
      custom_lenses;
      isp_meniscuses; % Defines all their positive and negative meniscus lenses
      
      % analysis parameters
      P.n_overlap_spots = ceil(3000/(2*P.mirror_spacing));
      P.beam_diameter = 0.4;
      
      P.Herriott_passes = 1000; % essentially unlimited. Stop with HR=0.
      P.stop_ICOS = 0; % makes mirror 3 black
      P.ICOS_passes_per_injection = 20;
      P.max_rays = 60;
      
      P.visible = false;
      P.edges_visible = true;
      P.visibility = [];
      P.view = [];
      P.plot_endpoints = 0;
      P.evaluate_endpoints = 0;
      P.skip.total_power = 0;
      P.skip.overlap = 0;
      P.skip.eccentricity = 0;
      P.skip.mean_angle = 0;
      P.skip.RIM_passes = 0;
      P.focus = 0;
      P.propagate = 1;
    end
    
    function M = P_model(P)
      T = P.T; % 250e-6; % transmittance
      n_ZnSe = P.optics_n; % 2.4361;
      n_air = 1;
      d = P.mirror_spacing;
      if P.focus == 1
        n_optics = 4 + length(P.Lenses);
      else
        n_optics = 3;
      end
      visibility = P.visible * ones(n_optics,1);
      if length(P.visibility) <= n_optics
        visibility(1:length(P.visibility)) = P.visibility;
      end
      M = opt_model(n_optics, P.max_rays);
      M.Spot_Size = P.beam_diameter;
      M.visible = P.visible;
      %M.Optic{1} = HRmirror('HM', P.Hr, P.HRC, CT, 0, P.HR, ...
      %  [-P.herriott_spacing 0 0], [1 0 0], n_air, n_air, P.visible);
      M.Optic{2} = HRmirror('M1', P.r1, P.R1, P.CT1, T, 1-T, [0 0 0], ...
        [1 0 0], n_ZnSe, n_air, P.visible && visibility(2));
      % Ignore rays transmitted back through ICOS mirror:
      % M.Optic{2}.Surface{2}.emission_threshold = 2*T^2;
      M.Optic{2}.max_passes = P.Herriott_passes;
      M.Optic{2}.edges_visible = P.edges_visible;
      if P.stop_ICOS
        M.Optic{3} = HRmirror('M2', P.r2, P.R2, P.CT2, 0, 0, [d 0 0], ...
          [-1 0 0], n_ZnSe, n_air, P.visible && visibility(3));
      else
        M.Optic{3} = HRmirror('M2', P.r2, P.R2, P.CT2, T, 1-T, [d 0 0], ...
          [-1 0 0], n_ZnSe, n_air, P.visible && visibility(3));
        M.Optic{3}.max_passes = P.ICOS_passes_per_injection;
      end
      M.Optic{3}.edges_visible = P.edges_visible;
      if P.focus > 0
        opt_n = 4;
        opt_X = d + P.CT2;
        for i = 1:length(P.Lenses)
          opt_X = opt_X + P.Lens_Space(i);
          L = P.LensTypes.(P.Lenses{i});
          if strcmp(L.type,'positive_meniscus')
            if L.EFL > 0
              D = [-1,0,0];
            else
              D = [1, 0, 0];
            end
            M.Optic{opt_n} = positive_meniscus(L.r, L.R1, L.R2, L.CT, L.EFL, ...
              sprintf('L%d', i), [opt_X,0,0], D, n_ZnSe, n_air, ...
              P.visible && visibility(opt_n));
          elseif strcmp(L.type,'negative_meniscus')
            M.Optic{opt_n} = negative_meniscus(L.r, L.R1, L.R2, L.CT, L.EFL, ...
              sprintf('L%d', i), [opt_X,0,0], [-1,0,0], n_ZnSe, n_air, ...
              P.visible && visibility(opt_n));
          end
          M.Optic{opt_n}.edges_visible = P.edges_visible;
          opt_X = opt_X + L.CT;
          opt_n = opt_n + 1;
        end
        opt_X = opt_X + P.detector_spacing;
        M.Optic{opt_n} = detector(P.D_l, [opt_X,P.D_dY,0], [-1,0,0], ...
          P.visible && visibility(opt_n), 15);
      end
      m = P.injection_scale;
      APincident = M.Optic{2}.O + [0, m*P.y0, m*P.z0];
      Dincident = [1, -m*P.dy, m*P.dz];
      Ap = APincident - P.herriott_spacing*Dincident;
      M.Optic{1} = Herriott_Mirror('HM', P.Hr, P.HRC, P.HCT, P.HR, Ap, ...
        P.beam_diameter, [-P.herriott_spacing 0 0], [1 0 0], ...
        P.visible && visibility(1));
      Pincident = M.Optic{2}.O + [0, m*P.y0 + P.beam_dy, m*P.z0 + P.beam_dz];
      Oincident = Pincident - (P.herriott_spacing+1)*Dincident;
      if P.propagate
        M.push_ray(opt_ray(Oincident, Dincident), 0, 1, 2);
      end
    end
    
    function Res = results_struct
      Res.total_power = [];
      Res.max_radius = [];
      Res.inside = [];
      Res.overlap = [];
      Res.eccentricity = [];
      Res.mean_angle = [];
      Res.n_rays = [];
      Res.max_rays = [];
      Res.RIM_passes = [];
    end
  end
end
