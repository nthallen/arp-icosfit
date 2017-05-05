classdef WhiteCell < opt_model_p
  % This model will include a White Cell
  % Optic:
  %  1: WhiteCell Mirror with two aperatures
  %  2: Alignment Mirror 1
  %  3: Alignment Mirror 2
  % The salient features include:
  %   Cell length
  %   Mirror radius of curvature (assumed the same for all mirrors,
  %     most likely matching the cell length)
  %   Alignment Mirror spacing (horizontal distance between centers)
  %   Alignment Mirror radius (less than 1/2 mirror spacing)
  %   Alignment Foci spacing (could be independently adjusted)
  %   WhiteCell Mirror radius
  %   Beam entry position, direction (direction could be derived so
  %     beam hits the center of the closer alignment mirror)
  %   Detector position, direction
  %------
  % WhiteCell Mirror geometry is a little tricky. It could be hacked by
  % using the Herriott Mirror and placing the detector in front of the
  % mirror at the exit position. Clearly a nicely implmented object with
  % the correct geometry would be great.
  
  
  properties
    % apparently all the important properties are maintained in the
    % opt_model_p and opt_model structures.
  end
  
  methods
    function PM = WhiteCell(P,varargin)
      PM = PM@opt_model_p(P,varargin{:});
    end

    function [optics,pass_out] = evaluate_passes(PM)
      rays = 1:PM.M.n_rays;
      zs = find([ PM.M.Rays(rays).n_inc ] == 0);
      pass = rays-interp1(zs,zs,rays,'previous','extrap')+1;
      % Now I want to calculate the spot size and power on each pass
      max_pass = max(pass);
      opt_n = zeros(max_pass,1);
      opt_n(2:2:max_pass) = 2;
      opt_n(4*PM.M.User.N4) = 1;
      opt_n(1:4:max_pass) = 3;
      opt_n(3:4:max_pass) = 4;
      optics = opt_n;
      if nargout > 1
        pass_out = pass;
      end
    end
   
    function Res = evaluate_endpoints(PM,P)
      Res = PM.results_struct();
      if P.evaluate_endpoints <= 0
        return;
      end
      
      Res = PM.results_struct;
      Res.Inside = 1;
      
      Res.max_rays = PM.M.max_rays;
      
      [opt_n,pass] = PM.evaluate_passes;
      max_pass = length(opt_n);
      
      if P.evaluate_endpoints > max_pass
        return;
      end
      
      rays = 1:PM.M.n_rays;
      n_opt = [PM.M.Rays(rays).n_opt];
      Ipower = 0;
      Opower = 0;
      sizes = 0;
      position = zeros(1,3);
      for i = P.evaluate_endpoints
        % fprintf(1,'Pass = %d\n', i);
        vf = find(pass == i);
        xyz = zeros(length(vf),3);
        for j=1:length(vf)
          vfi = vf(j);
          % fprintf(1, 'j = %d vf(j) = %d\n', j, vfi);
          ray = PM.M.Rays(vfi).ray;
          xyz(j,:) = ray.E;
          if opt_n(i) == n_opt(vfi) && ray.Inside
            Ipower = Ipower + ray.P;
          else
            Opower = Opower + ray.P;
          end
        end
        sizes = sqrt(var(xyz(:,2))+var(xyz(:,3)));
        position = mean(xyz) - PM.M.Optic{opt_n(i)}.O;
      end
      Ipower = Ipower/P.beam_samples;
      Opower = Opower/P.beam_samples;
      Res.Ipower = Ipower;
      Res.Opower = Opower;
      Res.spot_dia = sizes*2;
    end
 
  end
  
  methods (Static)
    function P = props(cell_len, N4)
      % The salient features include:
      %   Cell length
      %   Mirror radius of curvature (assumed the same for all mirrors,
      %     most likely matching the cell length)
      %   Alignment Mirror spacing (horizontal distance between centers)
      %   Alignment Mirror radius (less than 1/2 mirror spacing)
      %   Alignment Foci spacing (could be independently adjusted)
      %   WhiteCell Mirror radius
      %   Beam entry position, direction (direction could be derived so
      %     beam hits the center of the closer alignment mirror)
      %   Detector position, direction
      if nargin < 1
        P.Cell_Length = 25;
      else
        P.Cell_Length = cell_len;
      end
      if nargin < 2
        N4 = 5; % 20 passes default
      end
      
      P.N4 = N4;
      P.M0_r_max = 2.54;
      P.M0_r_min = 0.47*2.54; % 0.94" top to bottom
      P.M0_ap_dr = P.M0_r_min; % y length of notch
      P.M0_roc = P.Cell_Length;
      P.M0_CT = 0.5;
      P.M0_Ody = 0; % offset of center position
      P.M0_Odz = 0;
      P.M0_Ddy = 0; % offset of direction
      P.M0_Ddz = 0;
      
      P.y0 = P.M0_r_max - P.M0_ap_dr/2;
      P.z0 = P.M0_r_min/2;
      
      P.M1_r = 3*P.M0_r_max/8;
      P.M1_roc = P.M0_roc;
      P.M1_dy = +P.M1_r+.125*2.54; % position of M1
      P.M1_dz = 0;
      P.M1_fy = P.y0/(2*N4);
      P.M1_fz = 0;
      P.M1_CT = 0.5;
      P.M1_Ddz = 0;
      P.M12_Ddz = 0;
      
      P.M2_r = 3*P.M0_r_max/8;
      P.M2_roc = P.M1_roc;
      P.M2_dy = -P.M2_r-.125*2.54; % position of M2
      P.M2_dz = 0;
      P.M2_fy = -P.y0/(2*N4); % location of the focal point at M0
      P.M2_fz = 0;
      P.M2_CT = 0.5;
      P.M2_Ddz = 0;

      P.dy = (P.M1_dy - P.y0)/P.Cell_Length;
      P.dz = -P.z0/P.Cell_Length;
      P.beam_dy = 0;
      P.beam_dz = 0;
      P.dydy = 0;
      P.dzdz = 0;
      
      P.Det_l = 0.56;
      P.Det_x = -1.0;
      
      P.beam_diameter = 0.3;
      P.beam_samples = 100;
      P.beam_divergence = 1; % degrees
      P.max_rays = N4*5*P.beam_samples; % round up a bit
      P.T = 0; % No transmittance
      P.HR = 0.99; % Mirror reflectivity
      P.visible = false;
      P.edges_visible = true;
      P.visibility = [];
      P.view = [];
      P.plot_endpoints = 0;
      P.evaluate_endpoints = 0;
      P.propagate = 1;
      P.Tinterval = 60;
    end
    
    function M = P_model(P)
      n_optics = 4;
      visibility = P.visible * ones(n_optics,1);
      if length(P.visibility) <= n_optics
        visibility(1:length(P.visibility)) = P.visibility;
      end
      M = opt_model(n_optics, P.max_rays);
      M.User.N4 = P.N4;
      M.Spot_Size = P.beam_diameter;
      M.visible = P.visible;

      M0_Ap = [0,-P.y0,P.z0];
      M.Optic{2} = Herriott_Mirror('M0', P.M0_r_max, P.M0_roc, P.M0_CT, P.HR, ...
        M0_Ap, P.z0, [0,P.M0_Ody,P.M0_Odz], [1,P.M0_Ddy,P.M0_Ddz], ...
        P.visible && visibility(2));
      M.Optic{2}.alternate = true;

      MD = [-P.Cell_Length, P.M1_fy-P.M1_dy, P.M1_fz-P.M1_dz+P.M12_Ddz+P.M1_Ddz];
      MD = MD/sqrt(sum(MD.^2));
      M.Optic{3} = HRmirror('M1', P.M1_r, P.M1_roc, P.M1_CT, 0, P.HR, ...
        [P.Cell_Length, P.M1_dy, P.M1_dz], ...
        MD, 1, 1, ...
        P.visible && visibility(3));
      M.Optic{3}.alternate = true;
      M.Optic{3}.max_passes = 1000;
      
      MD = [-P.Cell_Length, P.M2_fy-P.M2_dy, P.M2_fz-P.M2_dz+P.M12_Ddz+P.M2_Ddz];
      MD = MD/sqrt(sum(MD.^2));
      M2_O = [P.Cell_Length,P.M2_dy,P.M2_dz];
      M.Optic{4} = HRmirror('M2', P.M2_r, P.M2_roc, P.M2_CT, 0, P.HR, ...
        M2_O, MD, 1, 1, ...
        P.visible && visibility(4));
      M.Optic{4}.max_passes = 1000;

      detD = M2_O - M0_Ap;
      detD = detD/sqrt(sum(detD.^2));
      detO = M0_Ap + P.Det_x * detD;
      M.Optic{1} = detector(P.Det_l,detO,detD,P.visible && visibility(1),45);
      
      if isfield(P,'rng_state')
        M.User.rng_state = P.rng_state;
        rng(M.User.rng_state);
      else
        M.User.rng_state = rng;
      end
      %M.Optic{4} = detector(P.Det_l, [P.Det_x,P.Det_y,P.Det_z], [1,0,0], ...
      %    P.visible && visibility(4), 15);
      if P.propagate
        if P.beam_samples > 1
          M = WhiteCell.WC_Sample(M,P);
        else
          Oincident = [0,P.y0+P.beam_dy,P.z0+P.beam_dz];
          Dincident = [1,P.dy+P.dydy,P.dz+P.dzdz];
          M.push_ray(opt_ray(Oincident, Dincident), 0, 1, 3);
        end
        rays = 1:M.n_rays;
        zs = find([ M.Rays(rays).n_inc ] == 0);
        M.User.pass = rays-interp1(zs,zs,rays,'previous','extrap')+1;
      end
    end
    
    function M = WC_Sample(M, P)
      Nsamples = P.beam_samples;
      X = rand(Nsamples,1);
      r = P.beam_diameter*sqrt(-log(X))/2;
      th = 2*pi*rand(Nsamples,1);
      Dy = cos(th);
      Dz = sin(th);
      dd = sind(P.beam_divergence)*r*2/P.beam_diameter;
      dy = r.*Dy;
      dz = r.*Dz;
      dydy = dd.*Dy;
      dzdz = dd.*Dz;
      TStart = tic;
      Treport = 0;
      for i = 1:Nsamples
        Oincident = [0, P.y0 + dy(i), P.z0 + dz(i)];
        Dincident = [1, P.dy + dydy(i), P.dz + dzdz(i)];
        M.push_ray(opt_ray(Oincident, Dincident), 0, 1, 3);
        M.propagate;
        TIter = toc(TStart);
        if TIter > Treport + P.Tinterval
          fprintf(1,'%.1f: Iteration: %d samples of %d\n', TIter, ...
            i, Nsamples);
          Treport = TIter;
        end
      end
    end
    
    function Res = results_struct
      Res.Ipower = [];
      Res.Opower = [];
      Res.spot_dia = [];
      Res.Inside = [];
      Res.max_rays = [];
    end
  end
end
