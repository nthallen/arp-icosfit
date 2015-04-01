classdef opt_model < handle
  properties
    Optic % linear array of optics
    n % index of refraction of air (may be function of nu, P)
    visible % whether or not to draw
    Rays
    n_rays
    Spot_Size % diameter in cm
  end
  
  properties
    max_rays
    ray_stack
    ray_stack_size
    ray_stack_idx
  end
  
  methods
    function M = opt_model(n_optics, max_rays)
      M.Optic = cell(n_optics,1);
      M.Rays = struct( ...
        'n_inc', 0, 'n_opt', 0, 'ray', cell(max_rays,1));
      M.visible = true;
      M.n_rays = 0;
      M.max_rays = max_rays;
      M.ray_stack_size = 10;
      M.ray_stack = struct( ...
        'n_inc', 0, ...
        'from_obj', 0, ...
        'to_obj', 0, ...
        'ray', cell(M.ray_stack_size,1));
      M.ray_stack_idx = 0;
      M.Spot_Size = 0.4; % 4 mm
    end
    
    function propagate(M)
      if M.visible
        clf;
      end
      for i=1:length(M.Optic)
        if isempty(M.Optic{i})
          error('MATLAB:HUARP:UndefinedOptic', 'Optic %d is undefined', i);
        end
        % M.Optic{i}.visible = M.visible;
        M.Optic{i}.draw;
      end
      if M.visible
        set(gca,'dataaspectratio',[1 1 1]);
      end
      while M.ray_stack_idx > 0 && M.n_rays < M.max_rays
        RS = M.pop_ray;
        Rincident = RS.ray;
        opt_n = RS.to_obj;
        if opt_n > 0 && opt_n <= length(M.Optic)
          [Rincident, Rreflect, ~, Rtransmit] = ...
            M.Optic{RS.to_obj}.propagate(Rincident);
          inc_n = M.save_ray(Rincident, RS.n_inc, opt_n);
          M.push_ray(Rreflect, inc_n, opt_n, RS.from_obj);
          if RS.from_obj < RS.to_obj
            next_opt = RS.to_obj + 1;
          else
            next_opt = RS.to_obj - 1;
          end
          M.push_ray(Rtransmit, inc_n, opt_n, next_opt);
        end
      end
    end
    
    function push_ray(M, R, n_inc, fr_obj, to_obj)
      if ~isempty(R) && to_obj > 0 && to_obj <= length(M.Optic)
        if M.ray_stack_idx >= M.ray_stack_size
          error('MATLAB:HUARP:RayStackOverflow', 'Ray Stack Overflow');
        end
        M.ray_stack_idx = M.ray_stack_idx + 1;
        M.ray_stack(M.ray_stack_idx).ray = R;
        M.ray_stack(M.ray_stack_idx).n_inc = n_inc;
        M.ray_stack(M.ray_stack_idx).from_obj = fr_obj;
        M.ray_stack(M.ray_stack_idx).to_obj = to_obj;
      end
    end
    
    function RS = pop_ray(M)
      if M.ray_stack_idx == 0
        error('MATLAB:HUARP:RayStackUnderflow', 'Ray Stack Underflow');
      end
      RS = M.ray_stack(M.ray_stack_idx);
      M.ray_stack_idx = M.ray_stack_idx - 1;
    end
    
    function n = save_ray(M, R, inc_n, opt_n)
      if M.n_rays >= M.max_rays
        error('MATLAB:HUARP:RayBufferOverflow', 'Ray Buffer Overflow');
      end
      M.n_rays = M.n_rays + 1;
      M.Rays(M.n_rays).n_inc = inc_n;
      M.Rays(M.n_rays).n_opt = opt_n;
      M.Rays(M.n_rays).ray = R;
      n = M.n_rays;
    end
    function [xyz, oxyz] = extract_endpoints(M, opt_n)
      % xyz = M.extract_endpoints(opt_n);
      % [xyz, oxyz] = extract_endpoints(M, opt_n);
      % Extract endpoint vectors and optionally their origin points.
      if opt_n < 1 || opt_n > length(M.Optic)
        error('MATLAB:HUARP:OpticOOR', 'Invalid Optic Number');
      end
      vf = find([M.Rays(1:M.n_rays).n_opt] == opt_n);
      xyz = zeros(length(vf),3);
      if nargout > 1
        oxyz = zeros(length(vf),3);
      end
      for i=1:length(vf)
        xyz(i,:) = M.Rays(vf(i)).ray.E;
        if nargout > 1
          oxyz(i,:) = M.Rays(vf(i)).ray.O;
        end
      end
    end
    
    function plot_endpoints(M, opt_n)
      xyz = M.extract_endpoints(opt_n);
      Opt = M.Optic{opt_n};
      if abs(dot([1 0 0],Opt.D)) < .99
        error('MATLAB:HUARP:RotatedOptic', ...
          'Rotated Optical elements not supported in plot_endpoints');
      end
      % figure;
      clf;
      p = M.Optic{opt_n}.Surface{1}.perimeter;
      plot(xyz(:,2), xyz(:,3), '.r', p(:,2), p(:,3), 'k');
      set(gca,'DataAspectRatio',[1 1 1 ]);
      xlabel('cm');
    end
    
    function d_angle = extract_angles(M, opt_n)
      [xyz, oxyz] = M.extract_endpoints(opt_n);
      ep_angle = atan2(xyz(:,2),xyz(:,3));
      o_angle = atan2(oxyz(:,2),oxyz(:,3));
      d_angle = ep_angle - o_angle;
      v = d_angle < 0;
      d_angle = d_angle + v*2*pi;
    end
    
    function Res = evaluate_endpoints(M, opt_n, n_pts)
      % Res = M.evaluate_endpoints(opt_n[, n_pts]);
      % Analysis includes:
      %   min_dist: Minimum spacing from first spot within n_pts
      %   min_dist_n: Spot number of minimum spacing
      %   area: area filled (assuming precession)
      %   mean_angle: in degrees
      %   n_rays
      xyz = M.extract_endpoints(opt_n);
      if nargin < 3
        n_pts = size(xyz,1);
      end
      
      if n_pts > 1
        d = xyz(2:n_pts,:) - ones(n_pts-1,1)*xyz(1,:);
        dd = sqrt(sum(d.^2,2));
        Res.min_dist = min(dd);
        Res.min_dist_n = find(dd == Res.min_dist,1);
        Res.close_n = { find(dd < M.Spot_Size) };
        Res.close = { dd(Res.close_n{1}) };
      else
        Res.min_dist = NaN;
        Res.min_dist_n = NaN;
        Res.close_n = NaN;
        Res.close = NaN;
      end
      rr = minmax(sqrt(sum(xyz(:,[2 3]).^2,2))');
      Res.area = pi*(rr(2)^2 - rr(1)^2);
      Res.mean_angle = mean(M.extract_angles(opt_n));
      if Res.mean_angle > pi
        Res.mean_angle = Res.mean_angle - 2*pi;
      end
      Res.mean_angle = rad2deg(Res.mean_angle);
      Res.n_rays = size(xyz,1);
      Res.max_rays = M.max_rays;
    end
  end
  
  methods (Static)
    function save_figure(f, fname)
      if strcmpi(fname(end-3:end),'.png')
        print(f, '-dpng', fname);
      elseif strcmpi(fname(end-3:end), '.fig')
        savefig(f, fname);
      else
        error('MATLAB:HUARP:UnkownFiletype', 'Unknown File type in save_figure');
      end
      delete(f);
    end
    
    function fld = extract_result(Res, fldname)
      fld = reshape([Res.(fldname)],size(Res,1),size(Res,2));
    end
    
    function Res = results_struct
        Res.min_dist = [];
        Res.min_dist_n = [];
        Res.close_n = [];
        Res.close = [];
        Res.area = [];
        Res.mean_angle = [];
        Res.n_rays = [];
        Res.max_rays = [];
    end
  end
end