classdef optic
  properties
    name % for error messages
    O % Origin
    D % Reference Direction vector
    Surface % two surface objects, perhaps cell array
    n_int % index of refraction of the optic
    n_ext % index of refraction of the air
    max_passes % limits propagation
    visible % boolean
%         If Ray direction dot optic direction < 0
%         we can assume the ray will hit Surface(1) first
%         Hence optic.Direction can be considered the
%         normal coming out of the optic at Surface(1).
  end
  methods
    function opt = optic(nm, O_in, D_in, ni_in, ne_in, vis)
      opt.name = nm;
      opt.O = O_in;
      opt.D = D_in;
      opt.n_int = ni_in;
      opt.n_ext = ne_in;
      opt.max_passes = 0; % no limit
      opt.visible = vis;
      opt.Surface = cell(2,1);
    end
    
    function [Rincident, Rreflect, Rinternal, Rtransmit] = propagate(opt, Rincident)
      opt.Surface{1}.visible = opt.visible;
      opt.Surface{2}.visible = opt.visible;
      if dot(Rincident.D, opt.D) < 0
        order = [1 2];
      else
        order = [2 1];
      end
      [Rincident,Rref1,Rtra1] = opt.Surface{order(1)}.propagate(Rincident);
      [Rtra1,Rref2,Rtransmit] = opt.Surface{order(2)}.propagate(Rtra1);
      [Rref2,Rref2ref1,Rref2tra1] = opt.Surface{order(1)}.propagate(Rref2);
      Rinternal = [Rtra1, Rref2];
      Rreflect = [Rref1; Rref2tra1];
      if ~isempty(Rref2ref1)
        warning('MATLAB:HUARP:MultIntRef', ...
          'Multiple internal reflections in optic %s', opt.name);
      end
    end
    
    function draw(opt)
      if (opt.visible)
        opt.Surface{1}.draw;
        hold on;
        opt.Surface{2}.draw;
        p1 = opt.Surface{1}.perimeter;
        p2 = flipud(opt.Surface{2}.perimeter);
        if size(p1,1)>size(p2,1)
          p1 = p1(1:size(p2,1),:);
        end
        pm = (p1+p2)/2;
        X = [p1(:,1) pm(:,1) p2(:,1)];
        Y = [p1(:,2) pm(:,2) p2(:,2)];
        Z = [p1(:,3) pm(:,3) p2(:,3)];
        surfl(X,Y,Z);
        shading flat;
      end
    end
  end
end