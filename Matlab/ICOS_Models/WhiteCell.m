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

   
    function Res = evaluate_endpoints(PM,P)
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
% 
%       if ~(P.skip.overlap && P.skip.eccentricity)
%         xyz = PM.M.extract_endpoints(opt_n);
%         Res.n_rays = size(xyz,1);
%       end
%       
       Res = PM.results_struct;
       Res.inside = PM.M.Inside;
% 
%       % Total power calculation
%       if ~P.skip.total_power
%         vf = find([PM.M.Rays(1:PM.M.n_rays).n_opt] == opt_n);
%         Res.total_power = 0;
%         for i = 1:length(vf)
%           ray = PM.M.Rays(vf(i)).ray;
%           if ray.Inside
%             Res.total_power = Res.total_power + ray.P;
%           end
%         end
%       end
% 
%       % Overlap calculation
%       % Now ignores opt_n and always analyzes for optic #2,
%       % the first ICOS mirror.
%       if ~P.skip.overlap
%         n_opt = [PM.M.Rays(1:PM.M.n_rays).n_opt];
%         v = n_opt == 2; % was opt_n
%         vi = find(v);
%         RIM_pass = zeros(length(vi),1);
%         for i=1:length(vi)
%           RIM_pass(i) = PM.M.Rays(vi(i)).ray.pass(1);
%         end
%         RIM_passes = unique(RIM_pass);
%         overlap = zeros(size(RIM_passes));
%         for j=1:length(RIM_passes)
%           ji = find(RIM_pass == RIM_passes(j));
%           n_pts = length(ji);
%           if n_pts > P.n_overlap_spots
%             n_pts = P.n_overlap_spots;
%             ji = ji(1:n_pts);
%           end
%           if n_pts > 1
%             P0 = PM.M.Rays(vi(ji(1))).ray.E(1,2:3);
%             for i=2:n_pts
%               d = PM.M.Rays(vi(ji(i))).ray.E(1,2:3) - P0;
%               d = sqrt(sum(d.^2,2));
%               overlap(j) = overlap(j) + ...
%                 sum(max(0,P.beam_diameter-d))/P.beam_diameter;
%             end
%           end
%         end
%         Res.overlap = mean(overlap);
%       end
%       
%       if ~P.skip.mean_angle
%         Res.mean_angle = mean(PM.M.extract_angles(opt_n));
%         if Res.mean_angle > pi
%           Res.mean_angle = Res.mean_angle - 2*pi;
%         end
%         Res.mean_angle = rad2deg(Res.mean_angle);
%       end
      
      Res.max_rays = PM.M.max_rays;
    end
 
  end
  
  methods (Static)
    function P = props
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
      P.Cell_Length = 25;
      P.M0_r_max = 2.54;
      P.M0_r_min = 0.75;
      P.M0_ap_dr = 0.75;
      P.M0_roc = P.Cell_Length;
      P.M0_CT = 0.5;
      P.M1_r = P.M0_r_max/2;
      P.M1_roc = P.M0_roc;
      P.M1_dy = +P.M1_r;
      P.M1_fy = 0.1;
      P.M1_CT = 0.5;
      P.M2_r = P.M0_r_max/2;
      P.M2_roc = P.M1_roc;
      P.M2_dy = -P.M2_r;
      P.M2_fy = -0.1;
      P.M2_CT = 0.5;
      P.Beam_y = P.M0_r_max - P.M0_ap_dr/2;
      P.Beam_z = P.M0_r_min/2;
      P.Beam_dy = (P.M1_dy - P.Beam_y)/P.Cell_Length;
      P.Beam_dz = -P.Beam_z/P.Cell_Length;
      P.Det_l = 0.1;
      P.Det_x = 0;
      P.Det_y = -P.Beam_y;
      P.Det_z = P.Beam_z;
      
      P.beam_diameter = 0.3;
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
      n_optics = 4;
      visibility = P.visible * ones(n_optics,1);
      if length(P.visibility) <= n_optics
        visibility(1:length(P.visibility)) = P.visibility;
      end
      M = opt_model(n_optics, P.max_rays);
      M.Spot_Size = P.beam_diameter;
      M.visible = P.visible;
      Ap = [0, P.Beam_y, P.Beam_z];
      M.Optic{1} = Herriott_Mirror('M0', P.M0_r_max, P.M0_roc, P.M0_CT, 1, ...
        Ap, 1.1 * P.beam_diameter/2, [0,0,0], [1,0,0], ...
        P.visible && visibility(1));
      MD = [-P.Cell_Length, P.M1_fy-P.M1_dy, 0];
      MD = MD/sqrt(sum(MD.^2));
      M.Optic{2} = HRmirror('M1', P.M1_r, P.M1_roc, P.M1_CT, 0, 1, ...
        [P.Cell_Length, P.M1_dy, 0], ...
        MD, 1, 1, ...
        P.visible && visibility(2));
      MD = [-P.Cell_Length, P.M2_fy-P.M2_dy, 0];
      MD = MD/sqrt(sum(MD.^2));
      M.Optic{3} = HRmirror('M2', P.M2_r, P.M2_roc, P.M2_CT, 0, 1, ...
        [P.Cell_Length, P.M2_dy, 0], ...
        MD, 1, 1, ...
        P.visible && visibility(3));
      M.Optic{4} = detector(P.Det_l, [P.Det_x,P.Det_y,P.Det_z], [1,0,0], ...
          P.visible && visibility(4), 15);
      if P.propagate
        Oincident = [0,P.Beam_y,P.Beam_z];
        Dincident = [1,P.Beam_dy,P.Beam_dz];
        M.push_ray(opt_ray(Oincident, Dincident), 0, 1, 2);
      end
    end
    
    function Res = results_struct
      Res.total_power = [];
      Res.max_radius = [];
      Res.inside = [];
      Res.overlap = [];
      Res.mean_angle = [];
      Res.n_rays = [];
      Res.max_rays = [];
    end
  end
end
