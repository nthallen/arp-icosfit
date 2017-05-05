classdef opt_model < handle
  properties
    Optic % linear array of optics
    n % index of refraction of air (may be function of nu, P)
    visible % whether or not to draw
    Rays
    n_rays
    Spot_Size % diameter in cm. Used for overlap estimation
    n_points_overlap
    Inside
    User
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
      M.n_points_overlap = 30;
      M.Inside = true;
    end
    
    function propagate(M)
      if M.visible && M.n_rays == 0
        clf;
      end
      for i=1:length(M.Optic)
        if isempty(M.Optic{i})
          error('MATLAB:HUARP:UndefinedOptic', 'Optic %d is undefined', i);
        end
        % M.Optic{i}.visible = M.visible;
        if M.visible && M.n_rays == 0
          M.Optic{i}.draw;
        end
      end
      if M.visible
        set(gca,'dataaspectratio',[1 1 1]);
      end
      inc_n = 0; % need a default
      while M.ray_stack_idx > 0 && M.n_rays < M.max_rays
        RS = M.pop_ray;
        % Rincident = RS.ray;
        opt_n = RS.to_obj;
        if opt_n > 0 && opt_n <= length(M.Optic)
          [Rincident, Rreflect, ~, Rtransmit] = ...
            M.Optic{opt_n}.propagate(RS.ray);
          if RS.from_obj < opt_n
            next_opt = opt_n + 1;
            prev_opt = opt_n - 1;
          else
            next_opt = opt_n - 1;
            prev_opt = opt_n + 1;
          end
          if M.Optic{opt_n}.alternate && isempty(Rincident)
            M.push_ray(RS.ray, inc_n, RS.from_obj, next_opt);
          else
            inc_n = M.save_ray(Rincident, RS.n_inc, opt_n);
            max_passes = M.Optic{opt_n}.max_passes;
            if RS.from_obj < opt_n % Forward transmission
              if max_passes > 0
                if ~isempty(Rreflect)
                  Rreflect.pass(end) = Rreflect.pass(end) + 1;
                  if Rreflect.pass(end) > max_passes
                    Rreflect = [];
                  end
                end
                if ~isempty(Rtransmit)
                  Rtransmit.pass(end+1) = 1;
                end
              end
            else % backwards transmission
              Rtransmit = []; % suppress all backwards transmissions
            end
            M.push_ray(Rreflect, inc_n, opt_n, prev_opt);
            M.push_ray(Rtransmit, inc_n, opt_n, next_opt);
          end
        end
      end
    end
    
    function redraw(M, ray_select)
      % M.redraw();
      % M.redraw(ray_select);
      %   ray_select is a row vector identifying which rays to draw
      if M.visible
        clf;
        for i=1:length(M.Optic)
          if isempty(M.Optic{i})
            error('MATLAB:HUARP:UndefinedOptic', 'Optic %d is undefined', i);
          end
          M.Optic{i}.draw;
        end
        set(gca,'dataaspectratio',[1 1 1]);
        % drawnow;
        if nargin < 2
          ray_select = 1:M.n_rays;
        end
        for i=ray_select
          M.Rays(i).ray.draw;
          % drawnow;
        end
      end
    end
    
    function set_visible(M, vis, opt_n)
      if nargin < 3
        opt_n = 1:length(M.Optic);
      end
      if iscolumn(opt_n)
        opt_n = opt_n';
      end
      if iscolumn(vis)
        vis = vis';
      end
      if length(opt_n) > 1 && length(vis) == 1
        vis = vis*ones(size(opt_n));
      end
      if ~isrow(opt_n) || ~isrow(vis) || length(opt_n) ~= length(vis)
        error('MATLAB:HUARP:argument mismatch', ...
          'vis and opt_n not the same length');
      end
      for i=1:length(opt_n)
        M.Optic{opt_n(i)}.visible = vis(i);
      end
      mvis = false;
      for i=1:length(M.Optic)
        if M.Optic{i}.visible
          mvis = true;
        end
      end
      M.visible = mvis;
    end
    
    function push_ray(M, R, n_inc, fr_obj, to_obj)
      % M.push_ray(R, n_inc, fr_obj, to_obj)
      % R: ray
      % n_inc: the index of the preceding ray in Rays
      % fr_obj: the source optic index
      % to_obj: the destination optic index
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
      if R.Inside == 0 % truly outside. Does not count exiting through
        % herriott aperature
        M.Inside = false;
      end
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
      drawnow; shg;
    end
    
    function plot_endpoints_skew(M, opt_n)
      Opt = M.Optic{opt_n};
      if abs(dot([1 0 0],Opt.D)) < .99
        error('MATLAB:HUARP:RotatedOptic', ...
          'Rotated Optical elements not supported in plot_endpoints');
      end
      [xyz, dxyz] = M.extract_endpoints_skew(opt_n);
      Y = [ xyz(:,2), xyz(:,2)+dxyz(:,2)];
      Z = [ xyz(:,3), xyz(:,3)+dxyz(:,3)];
      % figure;
      clf;
      p = M.Optic{opt_n}.Surface{1}.perimeter;
      plot(xyz(:,2), xyz(:,3), '.r', p(:,2), p(:,3), 'k', Y', Z', 'b');
      set(gca,'DataAspectRatio',[1 1 1 ]);
      xlabel('cm');
    end
    
    function [xyz, dxyz, divergence, skew] = extract_endpoints_skew(M, opt_n)
      % xyz in the point of incidence on the optic surface
      % dxyz is a direction vector normalized so the x component is 1
      [xyz, oxyz] = M.extract_endpoints(opt_n);
      dxyz = xyz - oxyz;
      % ldxyz = sqrt(sum(dxyz.^2,2));
      % dxyz = diag(1./ldxyz)*dxyz;
      dxyz = diag(1./dxyz(:,1)) * dxyz;
      
      if nargout >= 3
        yz = xyz(:,[2 3]);
        dyz = dxyz(:,[2 3]);
        r = sqrt(sum(yz.^2,2));
        ryz = diag(1./r) * yz; % unit vector in radial direction
        divergence = sum(ryz .* dyz,2);
      end
      if nargout >= 4
        col = ones(size(ryz,1),1);
        rsk = cross(col*[1 0 0], [0*col,ryz]);
        rsk = rsk(:,[2 3]);
        skew = sum(rsk.*dyz,2);
      end
    end
    
    function [oxyz, r, div, skew] = extract_origin_skew(M, opt_n)
      [xyz,oxyz] = M.extract_endpoints(opt_n);
      dxyz = xyz - oxyz;
      dxyz = diag(1./dxyz(:,1)) * dxyz;
      yz = oxyz(:,[2 3]); % position at origin, not dest.
      dyz = dxyz(:,[2 3]);
      r = sqrt(sum(yz.^2,2));
      ryz = diag(1./r) * yz; % unit vector in radial direction
      div = sum(ryz .* dyz,2);
      col = ones(size(ryz,1),1);
      rsk = cross(col*[1 0 0], [0*col,ryz]);
      rsk = rsk(:,[2 3]);
      skew = sum(rsk.*dyz,2);
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
      % This function is ripe for override to eliminate unnecesary
      % computation and to provide customized analysis.
      %
      % If you create new metrics and/or eliminate old ones, you
      % must also override the results_struct static method. Note
      % that these methods can also be overridden in subclasses of
      % opt_model_p.
      %
      % Current Analysis includes:
      %   min_dist: Minimum spacing from first spot within n_pts
      %   min_dist_n: Spot number of minimum spacing
      %   area: area filled (assuming precession)
      %   mean_angle: in degrees
      %   n_rays
      Res = M.results_struct;
      
      [xyz, oxyz] = M.extract_endpoints(opt_n);
      if nargin < 3
        n_pts = size(xyz,1);
      end

      % Overlap calculation
      Res.overlap = 0;
      if n_pts > 1
        if n_pts < M.n_points_overlap
          npo = n_pts;
        else
          npo = M.n_points_overlap;
        end
        for i=1:npo
          d = xyz(i+1:npo,:) - ones(npo-i,1)*xyz(i,:);
          d = sqrt(sum(d.^2,2));
          Res.overlap = Res.overlap + ...
            sum(max(0,M.Spot_Size-d))/M.Spot_Size;
        end
      end
      
      % inside calculation
      Res.inside = M.Inside;
      
      % Total power calculation
      vf = find([M.Rays(1:M.n_rays).n_opt] == opt_n);
      Res.total_power = 0;
      for i = 1:length(vf)
        Res.total_power = Res.total_power + M.Rays(i).ray.P;
      end
      
      % min_dist, min_dist_n, close and close_n
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
      Res.max_radius = rr(2);
      Res.eccentricity = sqrt(1 - (rr(1)^2 / rr(2)^2));
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
      Res.overlap = [];
      Res.inside = [];
      Res.n_rays = [];
      Res.max_rays = [];
      Res.max_radius = [];
      Res.eccentricity = [];
      Res.mean_angle = [];
      Res.total_power = [];

      Res.min_dist = [];
      Res.min_dist_n = [];
      Res.close_n = [];
      Res.close = [];
    end
  end
end
