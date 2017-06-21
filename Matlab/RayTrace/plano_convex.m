classdef plano_convex < optic
  properties
    r % radius
    R1 % convex surface radius of curvature
    CT % Center thickness
    EFL % Effective Focal Length?
  end
  
  methods
    function PM = plano_convex(r, R1, CT, EFL, ...
          nm, O_in, D_in, ni_in, ne_in, vis)
        % opt = double_convex(r, R1, CT, EFL, ...
        %   nm, O_in, D_in, ni_in, ne_in, vis);
        % D is normal out of the convex surface R1
      PM = PM@optic(nm, O_in, D_in, ni_in, ne_in, vis);
      PM.r = r;
      PM.CT = CT;
      PM.EFL = EFL; % 2*r;
      PM.R1 = R1; % M*PM.EFL; %
      PM.Surface{1} = spherical_segment(r, -PM.R1, PM.O, PM.D, 1, 0, PM.n_int, PM.n_ext);
      PM.Surface{2} = circular_surface(r, PM.O-PM.CT*PM.D, -PM.D, 1, 0, PM.n_int, PM.n_ext);
    end
  end
end
