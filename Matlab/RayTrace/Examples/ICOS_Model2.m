classdef ICOS_Model2 < opt_model_p
  properties
  end
  
  methods
    function PM = ICOS_Model2(P,varargin)
      PM = PM@opt_model_p(P,varargin{:});
    end
  end
  
  methods (Static)
    function P = props
      P.vis = false;
      P.d = 50; % mirror spacing
      P.r = 1.5*2.54; % mirror radius
      P.RC = 4.676*P.d;
      P.injection_scale = 1;
      P.max_rays = 60;
      P.y0 = 2.9;
      P.dy = .04;
      P.dz = .03;
    end
    
    function M = P_model(P)
      T = 250e-6; % transmittance
      CT = 0.2; % mirror center thickness
      n_ZnSe = 2.4361;
      n_air = 1;
      M = opt_model(2, P.max_rays);
      M.visible = P.vis;
      M.Optic{1} = HRmirror('M1', P.r, P.RC, CT, T, 1-T, [0 0 0], [1 0 0], n_ZnSe, n_air, P.vis);
      M.Optic{2} = HRmirror('M2', P.r, P.RC, CT, T, 1-T, [P.d 0 0], [-1 0 0], n_ZnSe, n_air, P.vis);
      D = [1, -P.dy * P.injection_scale, P.dz * P.injection_scale];
      P1 = [0, P.y0*P.injection_scale, 0];
      dX = -2;
      P0 = P1 + dX*D;
      M.push_ray(opt_ray(P0, D), 0, 0, 1);
    end
  end
end
