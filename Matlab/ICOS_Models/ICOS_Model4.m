classdef ICOS_Model4 < opt_model_p
  % This model will include a Herriot Mirror
  % Optic:
  %  1: Herriott Mirror
  %  2: HR mirror 1
  %  3: HR mirror 2
  %  4: Large focusing optic
  %  5: Small divergent focusing optic
  %  6: Small convergent focusing optic
  %  7: Detector
  properties
  end
  
  methods
    function PM = ICOS_Model4(P,varargin)
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
      vf = find([PM.M.Rays(1:PM.M.n_rays).n_opt] == opt_n);
      Res.total_power = 0;
      for i = 1:length(vf)
        ray = PM.M.Rays(vf(i)).ray;
        if ray.Inside
          Res.total_power = Res.total_power + ray.P;
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
      P.focus = 0;
      P.herriott_spacing = 10; % Before the first ICOS mirror
      P.HRC = 15*2.54; % Herriott radius of curvature
      P.Hr = 1.5*2.54; % Herriott radius
      P.HR = 0.98; % Herriott reflectivity
      P.ICOS_passes_per_injection = 20;
      P.stop_ICOS = 0;
      P.mirror_spacing = 50; % mirror spacing
      P.CT = 0.2; % mirror center thickness
      P.lens_spacing = 0.5;
      P.r = 1.5*2.54; % mirror radius
      P.lens2_spacing = 2.3 * P.r;
      P.RC = 4.676*P.mirror_spacing;
      P.max_rays = 60;
      P.injection_scale = 1; % Scales y0, dy and dz
      P.y0 = P.Hr-0.54; % Location of Herriott hole
      P.dy = .04;
      P.dz = .03;
      P.L1_R1 = 8.0122;
      P.L1_R2 = 29.8275;
      P.L1_CT = 0.9;
      P.L1_EFL = 7.62;
%       P.L2_r = 1.27/2;
%       P.L2_R1 = 1.079;
%       P.L2_R2 = 2.296;
%       P.L2_CT = 0.27;
%       P.L2_EFL = 1.3;
      % ZC-NM-25-25 1" negative meniscus d/f = 1
      P.LensTypes.ZC_NM_25_25.r = 2.54/2;
      P.LensTypes.ZC_NM_25_25.R1 = 2.366;
      P.LensTypes.ZC_NM_25_25.R2 = 6.427;
      P.LensTypes.ZC_NM_25_25.CT = 0.21;
      P.LensTypes.ZC_NM_25_25.EFL = 2.54;
      % ZC-PM-25-25 1" positive meniscus d/f = 1
      P.L3_r = 2.54/2;
      P.L3_R1 = 2.291;
      P.L3_R2 = 5.834;
      P.L3_CT = 0.44;
      P.L3_EFL = 2.54;
      % ZC-NM-25-25 1" negative meniscus d/f = 1
      P.L4_r = 2.54/2;
      P.L4_R1 = 2.366;
      P.L4_R2 = 6.427;
      P.L4_CT = 0.21;
      P.L4_EFL = 2.54;
      % ZC-PM-25-25 1" positive meniscus d/f = 1
      P.L5_r = 2.54/2;
      P.L5_R1 = 2.291;
      P.L5_R2 = 5.834;
      P.L5_CT = 0.44;
      P.L5_EFL = 2.54;
      
      P.detector_spacing = P.L2_EFL*2;
      P.L1_X = P.mirror_spacing + P.CT + P.lens_spacing;
      P.L2_X = P.L1_X + P.lens2_spacing;
      P.L3_X = P.L2_X + P.L2_CT + 0.2;
      P.L4_X = P.L3_X + 1.0; % arbitrary
      P.L5_dX = 0.2; % arbitrary
      P.D_X = P.L2_X + P.detector_spacing;
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
      LP = P.L1_X;
      L2P = P.L2_X;
      DP = P.D_X;
      if P.focus == 2
        n_optics = 9;
      elseif P.focus == 1
        n_optics = 7;
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
        M.Optic{4} = positive_meniscus(P.r, P.L1_R1, P.L1_R2, P.L1_CT, P.L1_EFL, ...
                'L1', [LP,0,0], [-1 0 0], n_ZnSe, n_air, P.visible && visibility(4));
        M.Optic{5} = negative_meniscus(P.L2_r, P.L2_R1, P.L2_R2, P.L2_CT, P.L2_EFL, ...
                'L2', [P.L2_X,0,0], [-1 0 0], n_ZnSe, n_air, P.visible && visibility(5));
        M.Optic{6} = positive_meniscus(P.L3_r, P.L3_R1, P.L3_R2, P.L3_CT, P.L3_EFL, ...
                'L3', [P.L3_X,0,0], [-1 0 0], n_ZnSe, n_air, P.visible && visibility(6));
        det_opt_n = 7;
        if P.focus == 2
          M.Optic{7} = negative_meniscus(P.L4_r, P.L4_R1, P.L4_R2, P.L4_CT, P.L4_EFL, ...
            'L4', [P.L4_X,0,0], [-1 0 0], n_ZnSe, n_air, P.visible && visibility(7));
          L5_X = P.L4_X + P.L4_CT + P.L5_dX;
          M.Optic{8} = positive_meniscus(P.L5_r, P.L5_R1, P.L5_R2, P.L5_CT, P.L5_EFL, ...
            'L5', [L5_X,0,0], [-1 0 0], n_ZnSe, n_air, P.visible && visibility(8));
          det_opt_n = 9;
        end
        M.Optic{det_opt_n} = detector(P.D_l, [DP,P.D_dY,0], [1,0,0],P.visible && visibility(det_opt_n));
      end
      m = P.injection_scale;
      Pincident = M.Optic{2}.O + [0, m*P.y0, 0];
      Dincident = [1, -m*P.dy, m*P.dz];
      Oincident = Pincident - (P.herriott_spacing+1)*Dincident;
      Ap = Pincident - P.herriott_spacing*Dincident;
      M.Optic{1} = Herriott_Mirror('HM', P.Hr, P.HRC, CT, P.HR, Ap, ...
        P.beam_diameter, [-P.herriott_spacing 0 0], [1 0 0], P.visible && visibility(1));
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
