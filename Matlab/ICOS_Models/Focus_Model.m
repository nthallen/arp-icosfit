classdef Focus_Model < opt_model_p
  % This model is just for looking at focus issues
  % Optic:
  %  1: Large focusing optic
  %  2: Detector
  properties
  end
  
  methods
    function PM = Focus_Model(P,varargin)
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
%       if n_pts > P.n_overlap_spots
%         n_pts = P.n_overlap_spots;
%       end
      
      Res = PM.results_struct;
      
%       Res.overlap = 0;
%       if n_pts > 1
%         for i=1:n_pts
%           d = xyz(i+1:n_pts,:) - ones(n_pts-i,1)*xyz(i,:);
%           d = sqrt(sum(d.^2,2));
%           Res.overlap = Res.overlap + ...
%             sum(max(0,P.beam_diameter-d))/P.beam_diameter;
%         end
%       end
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
    function [x, y, z] = identify_minimum(PM,metric)
      if ~isfield(PM.Results,metric)
        error('MATLAB:HUARP:NotAMetric', ...
          '"%s" is not a metric in this model', metric);
      end
      if isempty(PM.Results.x)
        error('MATLAB:HUARP:NoIndepententVariables', ...
          'There are no independent variables in this model');
      end
      z = nanmin(nanmin(PM.Results.(metric)));
      [i,j] = find(PM.Results.(metric) == z,1);
      x = PM.Results.(PM.Results.x)(i,j);
      if ~isempty(PM.Results.y)
        y = PM.Results.(PM.Results.y)(i,j);
      else
        y = [];
      end
    end
  end
  
  methods (Static)
    function P = props
      P.visible = true;
      P.plot_endpoints = 0;
      P.evaluate_endpoints = 0;
      P.max_rays = 60;
      P.injection_scale = 1; % Scales y0, dy and dz
      P.r = 1.5*2.54; % mirror radius
      P.y0 = 0; % Location of Herriott hole
      P.z0 = P.r;
      P.dy = 0;
      P.dz = 0;
      P.L1_R1 = 8.0122;
      P.L1_R2 = 29.8275;
      P.L1_CT = 0.9;
      P.L1_EFL = 7.62;
      P.lens2_spacing = P.L1_EFL;
      % This is the 1/2" focal length positive meniscus lens
      %       P.L2_r = 1.27/2;
      %       P.L2_R1 = 1.079;
      %       P.L2_R2 = 2.296;
      %       P.L2_CT = 0.27;
      %       P.L2_EFL = 1.3;
      % This is the 1" focal length positive meniscus lens
      P.L2_r = 2.54/2;
      P.L2_R1 = 2.366;
      P.L2_R2 = 6.427;
      P.L2_CT = 0.21;
      P.L2_EFL = 2.54;
      
      P.detector_spacing = P.L1_EFL;
      P.L1_X = 0;
      P.L2_X = P.L1_X + P.lens2_spacing;
      P.D_X = P.L1_X + P.detector_spacing;
      
      % analysis parameters
      % P.n_overlap_spots = ceil(3000/(2*P.mirror_spacing));
      % P.beam_diameter = 0.4;
    end
    
    function M = P_model(P)
      n_ZnSe = 2.4361;
      n_air = 1;
      n_optics = 2;
      M = opt_model(n_optics, P.max_rays);
      M.visible = P.visible;
      M.Optic{1} = positive_meniscus(P.r, P.L1_R1, P.L1_R2, P.L1_CT, P.L1_EFL, ...
              'L1', [P.L1_X,0,0], [-1 0 0], n_ZnSe, n_air, P.visible);
      % M.Optic{4} = positive_meniscus(P.L2_r, P.L2_R1, P.L2_R2, P.L2_CT, P.L2_EFL, ...
      %      'L2', [L2P,0,0], [1 0 0], n_ZnSe, n_air, P.visible);
      % M.Optic{5} = positive_meniscus(P.L2_r, P.L2_R1, P.L2_R2, P.L2_CT, P.L2_EFL, ...
      %         'L2', [L2P,0,0], [1 0 0], n_ZnSe, n_air, P.visible);
      M.Optic{2} = detector(.2, [P.D_X,0,0], [-1,0,0], P.visible);
      m = P.injection_scale;
      for m = linspace(0,0.95,6)
        L = 1+P.L1_CT;
        O = [-1, m*P.y0 + L*P.dy, m*P.z0 - L*P.dz];
        M.push_ray(opt_ray( O, [1, -P.dy, P.dz] ), 0, 0, 1);
      end
    end
    
    function Res = results_struct
      Res.max_radius = [];
      % Res.overlap = [];
      Res.eccentricity = [];
      Res.mean_angle = [];
      Res.n_rays = [];
      Res.max_rays = [];
    end
  end
end