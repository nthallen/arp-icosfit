classdef ICOS_Model3 < opt_model_p
  properties
  end
  
  methods
    function PM = ICOS_Model3(P,varargin)
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
  end
  
  methods (Static)
    function P = props
      P.visible = false;
      P.plot_endpoints = 0;
      P.mirror_spacing = 50; % mirror spacing
      P.CT = 0.2; % mirror center thickness
      P.lens_spacing = 0.5;
      P.r = 1.5*2.54; % mirror radius
      P.lens2_spacing = 2.3 * P.r;
      P.RC = 4.676*P.mirror_spacing;
      P.max_rays = 60;
      P.y0 = 2.9;
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
      P.L2_r = 2.54/2;
      P.L2_R1 = 2.366;
      P.L2_R2 = 6.427;
      P.L2_CT = 0.21;
      P.L2_EFL = 2.54;
      P.detector_spacing = P.L2_EFL*2;
      P.L1_X = P.mirror_spacing + P.CT + P.lens_spacing;
      P.L2_X = P.L1_X + P.lens2_spacing;
      P.D_X = P.L2_X + P.detector_spacing;
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
      M = opt_model(5, P.max_rays);
      M.visible = P.visible;
      M.Optic{1} = HRmirror('M1', P.r, P.RC, CT, T, 1-T, [0 0 0], [1 0 0], n_ZnSe, n_air, false);
      M.Optic{2} = HRmirror('M2', P.r, P.RC, CT, T, 1-T, [d 0 0], [-1 0 0], n_ZnSe, n_air, false);
      M.Optic{3} = positive_meniscus(P.r, P.L1_R1, P.L1_R2, P.L1_CT, P.L1_EFL, ...
              'L1', [LP,0,0], [-1 0 0], n_ZnSe, n_air, P.visible);
%     M.Optic{4} = positive_meniscus(P.L2_r, P.L2_R1, P.L2_R2, P.L2_CT, P.L2_EFL, ...
%             'L2', [L2P,0,0], [1 0 0], n_ZnSe, n_air, P.visible);
      M.Optic{4} = positive_meniscus(P.L2_r, P.L2_R1, P.L2_R2, P.L2_CT, P.L2_EFL, ...
              'L2', [L2P,0,0], [1 0 0], n_ZnSe, n_air, P.visible);
      M.Optic{5} = detector(6, [DP,0,0], [-1,0,0],P.visible);
      M.push_ray(opt_ray([-1, P.y0, 0], [1, -P.dy, P.dz]), 0, 0, 1);
    end
  end
end
