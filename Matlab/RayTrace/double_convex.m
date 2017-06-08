classdef double_convex < optic
  properties
    r % radius
    R1 % convex surface radius of curvature
    R2 % concave surface radius of curvature
    CT % Center thickness
    EFL % Effective Focal Length?
  end
  
  methods
    function PM = double_convex(r, R1, R2, CT, EFL, ...
          nm, O_in, D_in, ni_in, ne_in, vis)
        % D is normal out of R1
      PM = PM@optic(nm, O_in, D_in, ni_in, ne_in, vis);
      PM.r = r;
      PM.CT = CT;
      PM.EFL = EFL; % 2*r;
      PM.R1 = R1; % M*PM.EFL; %
      PM.R2 = R2; % (PM.n_int-1)*PM.EFL* ...
        % (1-(PM.CT*(PM.n_int-1)/(PM.M*PM.EFL*PM.n_int)))/...
        %  ((PM.n_int-1)/PM.M - 1);
      PM.Surface{1} = spherical_segment(r, -PM.R1, PM.O, PM.D, 1, 0, PM.n_int, PM.n_ext);
      PM.Surface{2} = spherical_segment(r, -PM.R2, PM.O-PM.CT*PM.D, -PM.D, 1, 0, PM.n_int, PM.n_ext);
    end
  end
end
