classdef Telescope < opt_model_p
  % Optical model of a telescope for investigating effects of divergence
  % Optic:
  %  1: Lens: ZC-PM-12-50
  %  2: Lens: ZC-PM-12-12
  %  3: Detector
  % The salient features include:
  %   Distance between elements
  
  
  properties
    % apparently all the important properties are maintained in the
    % opt_model_p and opt_model structures.
  end
  
  methods
    function PM = Telescope(P,varargin)
      PM = PM@opt_model_p(P,varargin{:});
    end

   
    function Res = evaluate_endpoints(PM,P,n_pts)
      Res = PM.results_struct();
      if P.evaluate_endpoints < 0
        return;
      elseif P.evaluate_endpoints > 0
        opt_n = P.evaluate_endpoints;
      elseif P.plot_endpoints > 0
        opt_n = P.plot_endpoints;
      else
        opt_n = length(PM.M.Optic);
      end
      Res = PM.results_struct;
      Res.inside = PM.M.Inside;
      Res.max_rays = PM.M.max_rays;
      [xyz, oxyz] = PM.M.extract_endpoints(opt_n);
      if nargin < 3
        n_pts = size(xyz,1);
      end
      
      % Total power calculation
      vf = find([PM.M.Rays(1:PM.M.n_rays).n_opt] == opt_n);
      Res.total_power = 0;
      for i = 1:length(vf)
        Res.total_power = Res.total_power + PM.M.Rays(i).ray.P;
      end
      
      % Max radius
      Res.max_radius = max(sqrt(sum(xyz(:,[2 3]).^2,2)));
      
      % Max divergence
      dxyz = xyz-oxyz;
      dxyz = diag(1./sqrt(sum(dxyz.^2,2)))*dxyz;
      ryz = diag(1./sqrt(sum(xyz(:,[2,3]).^2,2)))*(xyz*diag([0,1,1]));
      diverg = sum(dxyz.*ryz,2);
      Res.max_divergence = asind(nanmax(diverg));
    end
 
  end
  
  methods (Static)
    function P = props
      P.Lenses = { 'ZC_PM_12_50', 'ZC_PM_12_12' };
      P.L1_Space = 1; % Space before L1
      P.L2_Space = 6.35; % Space before L2
      P.optics_n = 2.4361;
      P.D_l = 0.2;
      P.D_Space = 1; % Space before Detector
      P.D_dY = 0;
      P.D_dZ = 0;
      P.beam_diameter = 0.4;
      % The following 7 parameters define an interface relevant to
      % ICOS_beam
      P.y0 = 0;
      P.z0 = 0;
      P.beam_dy = 0;
      P.beam_dz = 0;
      P.dy = 0;
      P.dz = 0;
      P.T = 1; % Transmittance, required by ICOS_beam
      % These three parameters invoke a different interface that provides
      % different visualization modes via opt_model and analysis via
      % opt_model_p. This interface is invoked when beam_n_th and beam_n_r
      % are non-zero.
      P.beam_divergence = 0;
      P.beam_n_th = 0;
      P.beam_n_r = 0;
      
      isp_meniscuses;
      
      P.max_rays = 60;
      P.visible = false;
      P.edges_visible = true;
      P.visibility = [];
      P.view = [];
      P.plot_endpoints = 0;
      P.evaluate_endpoints = 0;
      P.propagate = 1;
    end
    
    function M = P_model(P)
      n_optics = length(P.Lenses)+1;
      visibility = P.visible * ones(n_optics,1);
      if length(P.visibility) <= n_optics
        visibility(1:length(P.visibility)) = P.visibility;
      end
      lenses_n = P.optics_n;
      air_n = 1;
      % P.max_rays = 5;
      M = opt_model(n_optics, P.max_rays);
      M.Spot_Size = P.beam_diameter;
      M.visible = P.visible;
      opt_X = 0;
      opt_n = 1;
      Lens_Space = [P.L2_Space,P.D_Space];
      for i = 1:length(P.Lenses)
        L = P.LensTypes.(P.Lenses{i});
        if strcmp(L.type,'positive_meniscus')
          if L.EFL > 0
            D = [-1,0,0];
          else
            D = [1, 0, 0];
          end
          M.Optic{opt_n} = positive_meniscus(L.r, L.R1, L.R2, L.CT, L.EFL, ...
            sprintf('L%d', i), [opt_X,0,0], D, lenses_n, air_n, ...
            P.visible && visibility(opt_n));
        elseif strcmp(L.type,'negative_meniscus')
          M.Optic{opt_n} = negative_meniscus(L.r, L.R1, L.R2, L.CT, L.EFL, ...
            sprintf('L%d', i), [opt_X,0,0], [-1,0,0], lenses_n, air_n, ...
            P.visible && visibility(opt_n));
        end
        M.Optic{opt_n}.edges_visible = P.edges_visible;
        opt_X = opt_X + L.CT;
        opt_X = opt_X + Lens_Space(i);
        opt_n = opt_n + 1;
      end
      M.Optic{opt_n} = detector(P.D_l, [opt_X,P.D_dY,P.D_dZ], [-1,0,0], ...
        P.visible && visibility(opt_n), 15);
      if P.propagate
        if P.beam_n_r > 0 && P.beam_n_th > 0
        else
          Oincident = [-P.L1_Space,P.y0+P.beam_dy,P.z0+P.beam_dz];
          Dincident = [1,P.dy+P.dydy,P.dz+P.dzdz];
          Eincident = Oincident + P.L1_Space*Dincident;
          M.push_ray(opt_ray(Oincident, Dincident), 0, 0, 1);
        end
      end
    end
    
    function Res = results_struct
      Res.total_power = [];
      Res.max_radius = [];
      Res.max_divergence = [];
      Res.inside = [];
      Res.overlap = [];
      Res.mean_angle = [];
      Res.n_rays = [];
      Res.max_rays = [];
    end
  end
end
