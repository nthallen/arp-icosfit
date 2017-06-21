classdef WhiteCell < opt_model_p
  % This model will include a White Cell
  % Optic:
  %  1: Detector
  %  2: WhiteCell Mirror with two aperatures
  %  3: Alignment Mirror 1
  %  4: Alignment Mirror 2
  % Optics 170606:
  %  1: Detector
  %  2*: Post-LED Lens
  %  3*: Entrance window (notch) through M0
  %  4*: M0 Mirror
  %  5*: Exit window (notch)through M0
  %  6*: M1 Alignment Mirror
  %  7: M2 Alignment Mirror
  % Optics 170621:
  %  1: Detector
  %  2*: Post-LED Lens
  %  Housing:
  %    3*: Entrance aperture .25" cirular window
  %    4*: Housing wall circular dump
  %    5*: Exit aperture .25" circular window
  %  M0:
  %    6*: Entrance window (notch) through M0
  %    7*: M0 Mirror
  %    8*: Exit window (notch)through M0
  %  9*: M1 Alignment Mirror
  %  10: M2 Alignment Mirror
  %  
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
  %------
  % White Cell Mirror (M0) is a spherical mirror 2" in diameter that is
  % truncated 0.4375" above and below the horizontal diameter. There are
  % two symetrical rectangular notches demarcated by the path starting
  % where the horizontal diameter reaches the mirror's edge, extending in
  % along the diameter 0.266" and then up 0.4375" to the top of the mirror.
  %--
  % For the purposes of 170606, I will ignore the top and bottom
  % truncation, and model the mirror as a simple spherical mirror
  % overlapped by two square windows
  
  
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
      opt_n(1) = 2;
      opt_n(2) = 3;
      opt_n(3) = 6;
      opt_n(5:2:max_pass) = 7;
      opt_n(4:4:max_pass) = 9;
      opt_n(6:4:max_pass) = 10;
      opt_n(3+4*PM.M.User.N4) = 8;
      opt_n(4+4*PM.M.User.N4) = 5;
      opt_n(5+4*PM.M.User.N4) = 1;
      %opt_n(2:2:max_pass) = 2;
      %opt_n(4*PM.M.User.N4) = 1;
      %opt_n(1:4:max_pass) = 3;
      %opt_n(3:4:max_pass) = 4;
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
      
      % Res = PM.results_struct;
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
      P.M0_r_min = 0.4375*2.54; % (was 0.47:0.94" top to bottom)
      P.M0_Nl = 0.266*2.54; % y length of notch
      P.M0_roc = P.Cell_Length;
      P.M0_CT = 0.5;
      P.M0_Ody = 0; % offset of center position
      P.M0_Odz = 0;
      P.M0_Ddy = 0; % offset of direction
      P.M0_Ddz = 0;
      
      P.LED_dx = 0.3;
      P.Lens_dx = 3; % Distance from housing edge
      % 6mm x 9mm FL
      P.Lens_r = 0.1;
      P.Lens_CT = 0.27;
      P.Lens_ROC = 0.78;
      P.Lens_EFL = 0.9;
      P.Lens_type = 'double_convex';
      
      P.EnAp_y = P.M0_r_max - P.M0_Nl/2;
      P.EnAp_z = P.M0_Nl/2;
      P.EnAp_r = 0.125*2.54; % When I implement that
      
      P.Entrance_r = .125*2.54;
      P.Housing_r = P.M0_r_max+1;
      P.Housing_CT = (0.44 + 0.125)*2.54;
      
      P.M1_r = 3*P.M0_r_max/8;
      P.M1_roc = P.M0_roc;
      P.M1_dy = +P.M1_r+.125*2.54; % position of M1
      P.M1_dz = 0;
      P.M1_fy = P.EnAp_y/(2*N4);
      P.M1_fz = 0;
      P.M1_CT = 0.5;
      P.M1_Ddz = 0;
      P.M12_Ddz = 0;
      
      P.M2_r = 3*P.M0_r_max/8;
      P.M2_roc = P.M1_roc;
      P.M2_dy = -P.M2_r-.125*2.54; % position of M2
      P.M2_dz = 0;
      P.M2_fy = -P.EnAp_y/(2*N4); % location of the focal point at M0
      P.M2_fz = 0;
      P.M2_CT = 0.5;
      P.M2_Ddz = 0;

      % y0,z0,dy and dz will be determined by the EnAp_y, Ap_z, and
      % M1_dy. beam_d[yz] and d[yz]d[yz] can be used to explore
      % bad alignment.
      % P.y0 = 0; % derive from P.Ap_[yz]
      % P.z0 = 0; % derive from P.Ap_[yz]
      % P.dy = (P.M1_dy - P.y0)/P.Cell_Length;
      % P.dz = -P.z0/P.Cell_Length;
      P.beam_dy = 0;
      P.beam_dz = 0;
      P.dydy = 0;
      P.dzdz = 0;
      
      P.Det_l = 0.56;
      P.Det_dx = 0.25*2.54;
      
      P.beam_diameter = 0.17*2.54; % (was 0.3;)
      P.beam_samples = 100;
      P.beam_divergence = 6.2; % degrees
      % P.max_rays = N4*8*P.beam_samples; % round up a bit
      % P.T = 0; % No transmittance
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
      n_optics = 10;
      visibility = P.visible * ones(n_optics,1);
      if length(P.visibility) <= n_optics
        visibility(1:length(P.visibility)) = P.visibility;
      end
      M = opt_model(n_optics, P.N4*8*P.beam_samples);
      M.User.N4 = P.N4;
      M.Spot_Size = P.beam_diameter;
      M.visible = P.visible;

      M0_InAp = [-P.M0_CT,P.EnAp_y,P.EnAp_z]; % The inlet center
      M0_OutAp = [-P.M0_CT,-P.EnAp_y,P.EnAp_z]; % The inlet center
      
      M.Optic{3} = circular_window(P.EnAp_r, P.Housing_CT, M0_InAp, ...
        [1, 0, 0], P.visible && visibility(3));
      M.Optic{3}.alternate = true;
      M.Optic{4} = circular_dump(P.Housing_r, P.Housing_CT, ...
        [-P.M0_CT,0,0], [1, 0, 0], P.visible && visibility(4));
      M.Optic{4}.alternate = true;
      M.Optic{5} = circular_window(P.EnAp_r, P.Housing_CT, M0_OutAp, ...
        [1, 0, 0], P.visible && visibility(5));
      M.Optic{5}.alternate = true;
      M.Optic{5}.allow_reverse_transmission = true;
      
      
      M.Optic{7} = HRmirror('M0', P.M0_r_max, P.M0_roc, P.M0_CT, 0, P.HR, ...
        [0, P.M0_Ody,P.M0_Odz], [1,P.M0_Ddy,P.M0_Ddz], ...
        1, 1, P.visible && visibility(7));
      M.Optic{7}.alternate = true;
      
      Nh = sqrt(P.M0_Nl*(2-P.M0_Nl)); % Notch height
      Nl = max(P.M0_Nl,Nh);
      Ny = P.M0_r_max - P.M0_Nl + Nl/2;
      Nz = Nl/2;
      M.Optic{6} = square_window(Nl, P.M0_CT,[0,Ny,Nz],[1,0,0], ...
        P.visible && visibility(6));
      M.Optic{6}.alternate = true;
      M.Optic{8} = square_window(Nl, P.M0_CT,[0,-Ny,Nz],[1,0,0], ...
        P.visible && visibility(8));
      M.Optic{8}.alternate = true;
      M.Optic{8}.allow_reverse_transmission = true;
      

      MD = [-P.Cell_Length, P.M1_fy-P.M1_dy, P.M1_fz-P.M1_dz+P.M12_Ddz+P.M1_Ddz];
      MD = MD/sqrt(sum(MD.^2));
      M1_O = [P.Cell_Length, P.M1_dy, P.M1_dz];
      M.Optic{9} = HRmirror('M1', P.M1_r, P.M1_roc, P.M1_CT, 0, P.HR, ...
        M1_O, MD, 1, 1, ...
        P.visible && visibility(9));
      M.Optic{9}.alternate = true;
      M.Optic{9}.max_passes = 1000;
      
      MD = [-P.Cell_Length, P.M2_fy-P.M2_dy, P.M2_fz-P.M2_dz+P.M12_Ddz+P.M2_Ddz];
      MD = MD/sqrt(sum(MD.^2));
      M2_O = [P.Cell_Length,P.M2_dy,P.M2_dz];
      M.Optic{10} = HRmirror('M2', P.M2_r, P.M2_roc, P.M2_CT, 0, P.HR, ...
        M2_O, MD, 1, 1, ...
        P.visible && visibility(10));
      M.Optic{10}.max_passes = 1000;
      
      Lens_x = -P.M0_CT-P.Housing_CT-P.Lens_dx-P.Lens_CT;
      % M0_InAp is the center of the inside end of the enterance aperture 
      % M0_InAp1 is the center of the enterance aperature tube
      M0_InAp1 = M0_InAp-[P.Housing_CT/2,0,0];
      Lens_D = M0_InAp1-M1_O;
      % Lens_O is the lens origin
      Lens_O = M1_O + Lens_D*(M1_O(1)-Lens_x)/(M1_O(1)+P.M0_CT+P.Housing_CT/2);
      LED_O = M1_O + Lens_D*(M1_O(1)-Lens_x+P.LED_dx)/(M1_O(1)+P.M0_CT+P.Housing_CT/2);
      Lens_D = Lens_D/sqrt(sum(Lens_D.^2));
% ni fused silica @ 255 nm = 1.51
      if strcmp(P.Lens_type,'double_convex')
        M.Optic{2} = double_convex(P.Lens_r, P.Lens_ROC, P.Lens_ROC, ...
          P.Lens_CT, P.Lens_EFL, 'L1', Lens_O, Lens_D, 1.51, 1, ...
          P.visible && visibility(2));
      elseif strcmp(P.Lens_type,'plano_convex')
        M.Optic{2} = plano_convex(P.Lens_r, P.Lens_ROC, ...
          P.Lens_CT, P.Lens_EFL, 'L1', Lens_O, Lens_D, 1.51, 1, ...
          P.visible && visibility(2));
      else
        error('Invalid lens_type');
      end
      M.Optic{2}.alternate = true;
      
      LED_D = Lens_D/Lens_D(1); % denormalized to provide dy, dz

      Det_x = -P.Housing_CT-P.M0_CT-0.25*2.54;
      detD = M2_O - M0_OutAp;
      detD = detD/sqrt(sum(detD.^2));
      detO = M0_OutAp + Det_x * detD;
      M.Optic{1} = detector(P.Det_l,detO,detD,P.visible && visibility(1),45);
      
      if isfield(P,'rng_state')
        M.User.rng_state = P.rng_state;
        rng(M.User.rng_state);
      else
        M.User.rng_state = rng;
      end
      if P.propagate
        if P.beam_samples > 1
          M = WhiteCell.WC_Sample(M,P, LED_O, LED_D);
        else
          M.push_ray(opt_ray(LED_O, LED_D + [0,P.dydy,P.dzdz]), 0, 1, 2);
          M.propagate();
        end
        rays = 1:M.n_rays;
        zs = find([ M.Rays(rays).n_inc ] == 0);
        if length(zs) > 1
          M.User.pass = rays-interp1(zs,zs,rays,'previous','extrap')+1;
        else
          M.User.pass = rays;
        end
      end
    end
    
    function M = WC_Sample(M, P, LED_O, LED_D)
      Nsamples = P.beam_samples;
      X = rand(Nsamples,1);
      r = P.beam_diameter*sqrt(-log(X))/2;
      th = 2*pi*rand(Nsamples,1);
      Dy = cos(th);
      Dz = sin(th);
      dd = sind(P.beam_divergence)*r*2/P.beam_diameter;
      dy = r.*Dy;
      dz = r.*Dz;
      dydy = dd.*Dy + P.dydy;
      dzdz = dd.*Dz + P.dzdz;
      TStart = tic;
      Treport = 0;
      for i = 1:Nsamples
        Oincident = LED_O + [0, dy(i), dz(i)];
        Dincident = LED_D + [0, dydy(i), dzdz(i)];
        M.push_ray(opt_ray(Oincident, Dincident), 0, 1, 2);
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
