classdef ICOS_Model5a < opt_model_p
  % This model will include a Herriot Mirror
  % Optic:
  %  1: Herriott Mirror
  %  2: HR mirror 1
  %  3: HR mirror 2
  %  *: Lenses!
  %  x: Detector
  % This model differs from ICOS_Model5 by how the injection ray is
  % positioned:
  %
  % In ICOS_Model5:
  %   The position on the back of the first ICOS mirror (optic 2)
  %   is (0, y0, 0) and the direction is (1, -m*dy, m*dz). The
  %   ray's origin is back-calculated, and the position of the
  %   hole in the herriott mirror is defined at the intersection
  %   with the herriott's plane. This makes is a little easier
  %   to explore the (dy,dz) space without changing the point of
  %   insertion by much.
  % ICOS_Model5a:
  %   The position of the hole in the Herriott mirror is fixed.
  %   The direction is the same. This is useful for creating targets.
  properties
  end
  
  methods
    function PM = ICOS_Model5a(P,varargin)
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
      elseif P.plot_r_angle
        [xyz,dxyz] = PM.M.extract_endpoints_skew(4+length(P.Lenses));
        ldxyz = sqrt(sum(dxyz.^2,2));
        dxyz = diag(1./ldxyz)*dxyz;
        angle = rad2deg(acos(dxyz(:,1)));
        r = sqrt(sum(xyz(:,[2 3]).^2,2));
        plot(r, angle,'*');
        PM.draw_iteration_title(P);
        drawn = true;
      end
    end
    
    function Res = evaluate_endpoints(PM,P)
      if P.evaluate_endpoints > 0
        opt_n = P.evaluate_endpoints;
      elseif P.plot_endpoints > 0
        opt_n = P.plot_endpoints;
      else
        opt_n = length(PM.M.Optic);
      end
      xyz = PM.M.extract_endpoints(opt_n);
      n_pts = size(xyz,1);
      if n_pts > P.n_overlap_spots
        n_pts = P.n_overlap_spots;
      end
      
      Res = PM.results_struct;
      Res.inside = PM.M.Inside;

      % Total power calculation
      tan_angle_sq = [];
      if isfield(P,'acceptance_angle') && ~isempty(P.acceptance_angle)
        tan_angle_sq = tan(deg2rad(P.acceptance_angle))^2;
      end
      vf = find([PM.M.Rays(1:PM.M.n_rays).n_opt] == opt_n);
      Res.total_power = 0;
      for i = 1:length(vf)
        ray = PM.M.Rays(vf(i)).ray;
        if ray.Inside
          angle_ok = 1;
          if ~isempty(tan_angle_sq)
            dyz = ray.D([2 3])/ray.D(1);
            if sum(dyz.^2) > tan_angle_sq
              angle_ok = 0;
            end
          end
          if angle_ok
            Res.total_power = Res.total_power + ray.P;
          end
        end
      end

      % Overlap calculation
      Res.overlap = 0;
      if n_pts > 1
        for i=1:n_pts
          d = xyz(i+1:n_pts,:) - ones(n_pts-i,1)*xyz(i,:);
          d = sqrt(sum(d.^2,2));
          Res.overlap = Res.overlap + ...
            sum(max(0,P.beam_diameter-d))/P.beam_diameter;
        end
      end
      yz = xyz(:,[2 3]);
      rr = minmax(sqrt(sum(yz.^2,2))');
      Res.max_radius = rr(2);
      Res.eccentricity = sqrt(1 - (rr(1)^2 / rr(2)^2));
      Res.mean_angle = mean(PM.M.extract_angles(opt_n));
      if Res.mean_angle > pi
        Res.mean_angle = Res.mean_angle - 2*pi;
      end
      Res.mean_angle = rad2deg(Res.mean_angle);
      
      Res.n_rays = size(xyz,1);
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
      z = invert * z;
      x = PM.Results.(PM.Results.x)(i,j);
      if ~isempty(PM.Results.y)
        y = PM.Results.(PM.Results.y)(i,j);
      else
        y = [];
      end
      if nargout > 3
        ri = i;
      end
      if nargout > 4
        ci = j;
      end
    end
  end
  
  methods (Static)
    function P = props
      P.visible = false;
      P.visibility = [];
      P.plot_endpoints = 0;
      P.evaluate_endpoints = 0;
      P.plot_r_angle = 0;
      P.focus = 0;
      P.acceptance_angle = [];
      P.herriott_spacing = 10; % Before the first ICOS mirror
      P.HRC = 15*2.54; % Herriott radius of curvature
      P.Hr = 1.5*2.54; % Herriott radius
      P.HR = 0.98; % Herriott reflectivity
      P.ICOS_passes_per_injection = 20;
      P.stop_ICOS = 0;
      P.mirror_spacing = 50; % mirror spacing
      P.CT = 0.2; % mirror center thickness
      P.r = 1.5*2.54; % mirror radius
      P.RC = 4.676*P.mirror_spacing;
      P.max_rays = 60;
      P.injection_scale = 1; % Scales y0, dy and dz
      P.y0 = P.Hr-0.54; % Location of Herriott hole
      P.dy = .04;
      P.dz = .03;
      % Custom 3" poitive meniscus
      P.LensTypes.Lens1.type = 'positive_meniscus';
      P.LensTypes.Lens1.r = 3*2.54/2;
      P.LensTypes.Lens1.R1 = 8.0122;
      P.LensTypes.Lens1.R2 = 29.8275;
      P.LensTypes.Lens1.CT = 0.9;
      P.LensTypes.Lens1.EFL = 7.62;
      % ZC-NM-25-25 1" negative meniscus d/f = 1
      P.LensTypes.ZC_NM_25_25.type = 'negative_meniscus';
      P.LensTypes.ZC_NM_25_25.r = 2.54/2;
      P.LensTypes.ZC_NM_25_25.R1 = 2.366;
      P.LensTypes.ZC_NM_25_25.R2 = 6.427;
      P.LensTypes.ZC_NM_25_25.CT = 0.21;
      P.LensTypes.ZC_NM_25_25.EFL = 2.54;
      % ZC-NM-12-12 1/2" negative meniscus d/f = 1
      P.LensTypes.ZC_NM_12_12.type = 'negative_meniscus';
      P.LensTypes.ZC_NM_12_12.r = 2.54/4;
      P.LensTypes.ZC_NM_12_12.R1 = 1.072;
      P.LensTypes.ZC_NM_12_12.R2 = 2.938;
      P.LensTypes.ZC_NM_12_12.CT = 0.22;
      P.LensTypes.ZC_NM_12_12.EFL = -1.30;
      % ZC-NM-25-100 1" negative meniscus f = -10cm
      P.LensTypes.ZC_NM_25_100.type = 'negative_meniscus';
      P.LensTypes.ZC_NM_25_100.r = 2.54/2;
      P.LensTypes.ZC_NM_25_100.R1 = 10.424;
      P.LensTypes.ZC_NM_25_100.R2 = 41.4;
      P.LensTypes.ZC_NM_25_100.CT = 0.3;
      P.LensTypes.ZC_NM_25_100.EFL = -10;
      
      % ZC-PM-25-25 1" positive meniscus d/f = 1
      P.LensTypes.ZC_PM_25_25.type = 'positive_meniscus';
      P.LensTypes.ZC_PM_25_25.r = 2.54/2;
      P.LensTypes.ZC_PM_25_25.R1 = 2.291;
      P.LensTypes.ZC_PM_25_25.R2 = 5.834;
      P.LensTypes.ZC_PM_25_25.CT = 0.44;
      P.LensTypes.ZC_PM_25_25.EFL = 2.54;
      % ZC-PM-12-12 1/2" positive meniscus d/f = 1
      P.LensTypes.ZC_PM_12_12.type = 'positive_meniscus';
      P.LensTypes.ZC_PM_12_12.r = 2.54/4;
      P.LensTypes.ZC_PM_12_12.R1 = 1.079;
      P.LensTypes.ZC_PM_12_12.R2 = 2.296;
      P.LensTypes.ZC_PM_12_12.CT = 0.27;
      P.LensTypes.ZC_PM_12_12.EFL = 1.30;
      % ZC-PM-25-38 1" positive meniscus d/f = 1
      P.LensTypes.ZC_PM_25_38.type = 'positive_meniscus';
      P.LensTypes.ZC_PM_25_38.r = 2.54/2;
      P.LensTypes.ZC_PM_25_38.R1 = 3.804;
      P.LensTypes.ZC_PM_25_38.R2 = 12.273;
      P.LensTypes.ZC_PM_25_38.CT = 0.36;
      P.LensTypes.ZC_PM_25_38.EFL = 3.81;
      % ZC_PM_25_50 1" positive meniscus d/f = 0.5
      P.LensTypes.ZC_PM_25_50.type = 'positive_meniscus';
      P.LensTypes.ZC_PM_25_50.r = 2.54/2;
      P.LensTypes.ZC_PM_25_50.R1 = 4.395;
      P.LensTypes.ZC_PM_25_50.R2 = 10.765;
      P.LensTypes.ZC_PM_25_50.CT = 0.32;
      P.LensTypes.ZC_PM_25_50.EFL = 5.08;
      % ZC_PM_25_63 1" positive meniscus
      P.LensTypes.ZC_PM_25_63.type = 'positive_meniscus';
      P.LensTypes.ZC_PM_25_63.r = 2.54/2;
      P.LensTypes.ZC_PM_25_63.R1 = 5.849;
      P.LensTypes.ZC_PM_25_63.R2 = 16.108;
      P.LensTypes.ZC_PM_25_63.CT = 0.29;
      P.LensTypes.ZC_PM_25_63.EFL = 6.35;
      % ZC_PM_25_76 1" positive meniscus
      P.LensTypes.ZC_PM_25_76.type = 'positive_meniscus';
      P.LensTypes.ZC_PM_25_76.r = 2.54/2;
      P.LensTypes.ZC_PM_25_76.R1 = 6.823;
      P.LensTypes.ZC_PM_25_76.R2 = 17.947;
      P.LensTypes.ZC_PM_25_76.CT = 0.28;
      P.LensTypes.ZC_PM_25_76.EFL = 7.62;
      % ZC_PM_25_100 1" positive meniscus
      P.LensTypes.ZC_PM_25_100.type = 'positive_meniscus';
      P.LensTypes.ZC_PM_25_100.r = 2.54/2;
      P.LensTypes.ZC_PM_25_100.R1 = 9.078;
      P.LensTypes.ZC_PM_25_100.R2 = 25.000;
      P.LensTypes.ZC_PM_25_100.CT = 0.26;
      P.LensTypes.ZC_PM_25_100.EFL = 10.1;
      % ZC-PM-12-12 1/2" positive meniscus d/f = 1
      P.LensTypes.ZC_PM_12_12.type = 'positive_meniscus';
      P.LensTypes.ZC_PM_12_12.r = 2.54/4;
      P.LensTypes.ZC_PM_12_12.R1 = 1.079;
      P.LensTypes.ZC_PM_12_12.R2 = 2.296;
      P.LensTypes.ZC_PM_12_12.CT = 0.27;
      P.LensTypes.ZC_PM_12_12.EFL = 1.30;
      % ZC-HS-3-3 3mm hemisphere
      P.LensTypes.ZC_HS_3_3.type = 'positive_meniscus';
      P.LensTypes.ZC_HS_3_3.r = 0.15;
      P.LensTypes.ZC_HS_3_3.R1 = 0.151;
      P.LensTypes.ZC_HS_3_3.R2 = 10;
      P.LensTypes.ZC_HS_3_3.CT = 0.3;
      P.LensTypes.ZC_HS_3_3.EFL = -1.30;

      P.Lenses = { 'Lens1' };
      P.Lens_Space = [0.2];
      
      P.detector_spacing = 3;
      P.D_l = .2; % detector edge size
      P.D_dY = 0; % detector horizontal displacement
      
      % analysis parameters
      P.n_overlap_spots = ceil(3000/(2*P.mirror_spacing));
      P.beam_diameter = 0.4;
    end
    
    function M = P_model(P)
      T = 250e-6; % transmittance
      n_ZnSe = 2.4361;
      n_air = 1;
      CT = P.CT; % mirror center thickness
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
      m = P.injection_scale;
      Ap = [-P.herriott_spacing, P.y0, 0]; % Herriott mirror aperature
      M.Optic{1} = Herriott_Mirror('HM', P.Hr, P.HRC, CT, P.HR, Ap, ...
        P.beam_diameter, [-P.herriott_spacing 0 0], [1 0 0], P.visible && visibility(1));
      M.Optic{2} = HRmirror('M1', P.r, P.RC, CT, T, 1-T, [0 0 0], [1 0 0], n_ZnSe, n_air, P.visible && visibility(2));
      % Ignore rays transmitted back through ICOS mirror:
      M.Optic{2}.Surface{2}.emission_threshold = 2*T^2;
      if P.stop_ICOS
        M.Optic{3} = HRmirror('M2', P.r, P.RC, CT, 0, 0, [d 0 0], [-1 0 0], n_ZnSe, n_air, P.visible && visibility(3));
      else
        M.Optic{3} = HRmirror('M2', P.r, P.RC, CT, T, 1-T, [d 0 0], [-1 0 0], n_ZnSe, n_air, P.visible && visibility(3));
        % Limit the number of passes in the ICOS cell:
        M.Optic{3}.Surface{1}.emission_threshold = T*(1-T)^P.ICOS_passes_per_injection;
      end
      if P.focus > 0
        opt_n = 4;
        opt_X = d + P.CT;
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
              sprintf('L%d', i), [opt_X,0,0], D, n_ZnSe, n_air, P.visible && visibility(opt_n));
          elseif strcmp(L.type,'negative_meniscus')
            M.Optic{opt_n} = negative_meniscus(L.r, L.R1, L.R2, L.CT, L.EFL, ...
              sprintf('L%d', i), [opt_X,0,0], [-1,0,0], n_ZnSe, n_air, P.visible && visibility(opt_n));
          end
          opt_X = opt_X + L.CT;
          opt_n = opt_n + 1;
        end
        opt_X = opt_X + P.detector_spacing;
        M.Optic{opt_n} = detector(P.D_l, [opt_X,P.D_dY,0],[1,0,0],P.visible && visibility(opt_n));
      end
      Dincident = [1, -m*P.dy, m*P.dz];
      Oincident = Ap - (CT+2)*Dincident;
      M.push_ray(opt_ray(Oincident, Dincident), 0, 1, 2);
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
    end
  end
end
