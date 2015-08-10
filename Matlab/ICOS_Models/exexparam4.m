function res2 = exexparam4(loadfile)
% res2 = exexparam4(loadfile)
% loadfile is the name of a .mat file containing the output of
% exexparam3 or similar functions defining ICOS/RIM
% configuration parameters. This function extends those parameters
% to pick optimal focusing parameters.
%
% The output structure contains all the information from the input
% structure but adds additional fields for the optimized parameters.
%%
% load('exexparam3save.mat');
%%
load(loadfile);
%%
n = 2.4361;
i = 0;
P = ICOS_Model6.props;
LTS = fields(P.LensTypes);
clear res2
results = 0;
%%

%%
while i < length(res)
  %%
  i = i+1;
  %%
  P = render_model(res(i), 'visibility', [0 0 0], ...
    'focus', 1, 'ICOS_passes_per_injection', 100, ...
    'max_rays', 3000);
%   P.R1 = res(i).R1;
%   P.R2 = res(i).R2;
%   P.r1 = 3*2.54/2;
%   P.r2 = res(i).D2*2.54/2;
%   P.mirror_spacing = res(i).L;
%   P.y0 = res(i).Rr2;
%   P.dy = -res(i).ORd2;
%   P.dz = res(i).ORs2;
%   % P.CT2 = res(i).CT2;
%   P.HRC = res(i).RR1;
%   P.Hr = res(i).RD1*2.54/2;
%   P.HCT = 0.4;
%   P.herriott_spacing = res(i).ORL;
% 
%   P.stop_ICOS = 0;
%   P.visible = 1;
%   P.visibility = [0 0 0];
%   P.focus = 1;
%   P.ICOS_passes_per_injection = 100;
%   P.max_rays = 3000;
%   P.injection_scale = 1;
  % PM = ICOS_Model6(P);
  %%
  d = res(i).d2*n;
  s = res(i).s2;
  r = res(i).r2;
  j = 0;
  %%
  Lens1 = '';
  Li = 1;
  if d >= 0
    % place a PM lens at 0.2 that has a large enough diameter
    % and will give us a negative divergence.
    %%
    x = 0.2;
    rx = sqrt((r+d.*x).^2 + s.^2.*x.^2);
    dx = (d.*r + (d.^2+s.^2).*x)./rx;
    sx = s.*r./rx;
    f = P.LensTypes.Lens1.EFL;
    dxp = dx - rx/f;
    if dxp < 0
      Lens1 = 'Lens1 @ 0.2 + ';
    else
      Lens1 = 'ERROR ';
    end
    r = rx;
    d = dxp;
    s = sx;
    P.Lenses = {'Lens1'};
    P.Lens_Space = [0.2];
    Li = 2;
    %% recalculate d and continue
    P.HR = 0;
    P.evaluate_endpoints = 5;
    P.visible = 0;
    PM = ICOS_Model6(P);
    [xyz,oxyz] = PM.M.extract_endpoints(5);
    dxyz = xyz-oxyz;
    dxyz = diag(1./dxyz(:,1)) * dxyz;
    %%
    ro = sqrt(sum(oxyz(:,2:3).^2,2));
    dyz = sum(dxyz(:,2:3).*oxyz(:,2:3),2)./ro;
    rx = mean(ro);
    d = mean(dyz);
    s = s.*r./rx;
    r = rx;
  end
  %%
  while j < length(LTS)
    %%
    j = j+1;
    %%
    P.Lenses{Li} = LTS{j};
    Lens = P.LensTypes.(LTS{j});
    if strcmp(Lens.type,'positive_meniscus')
      f = Lens.EFL;
    else
      f = -Lens.EFL;
    end
    det_acc_limit = 14.9; % degrees acceptance angle
    [x,theta,ds] = pick_lens_x(r,d,s,f,det_acc_limit,Lens.r);
    %%
    if isempty(x)
%       P.Lens_Space(Li) = 1;
%       PM = ICOS_Model6(P);
%       title(sprintf('%d: Failure', i));
%       pause;
    else
      %%
      P.Lens_Space(Li) = x;
      P.detector_spacing = ds;
      P.visible = 0;
      PM = ICOS_Model6(P);
      %%
      [xyz,~,div,skew] = PM.M.extract_endpoints_skew(4+Li);
      theta2 = mean(rad2deg(atan(sqrt(div.^2 + skew.^2))));
      if theta2 > det_acc_limit + .01
%         warning('Matlab:HUARP:Angle', ...
%           'Calculated angle %.2f exceeds detector acceptance angle %.2f', ...
%           theta2, det_acc_limit);
        ls0 = P.Lens_Space(Li);
        theta0 = theta2;
        delta_ds = 0.4;
        theta1 = theta0;
        while theta1 <= theta0 && theta1 >= det_acc_limit
          ls1 = ls0+delta_ds;
          P.Lens_Space(Li) = ls1;
          PM = ICOS_Model6(P);
          [xyz,~,div,skew] = PM.M.extract_endpoints_skew(4+Li);
          theta1 = mean(rad2deg(atan(sqrt(div.^2 + skew.^2))));
          delta_ds = delta_ds + 0.1;
        end
        %%
        theta2 = theta1;
        if theta2 < det_acc_limit
          while abs(theta2-det_acc_limit) > .01
            ls2 = mean([ls0, ls1]);
            P.Lens_Space(Li) = ls2;
            PM = ICOS_Model6(P);
            [xyz,~,div,skew] = PM.M.extract_endpoints_skew(4+Li);
            theta2 = mean(rad2deg(atan(sqrt(div.^2 + skew.^2))));
            if theta2 > det_acc_limit
              ls0 = ls2;
            else
              ls1 = ls2;
            end
          end
        end
      end
      if theta2 <= det_acc_limit + .01
        %% Now adjust detector spacing using div and skew
        div = mean(div);
        skew = mean(skew);
        rx = mean(sqrt(sum(xyz(:,2:3).^2,2)));
        d_ds = -rx*div/(div^2+skew^2);
        ds = ds + d_ds;
        P.detector_spacing = ds;
        %%
%         P.evaluate_endpoints = 2;
%         P.skip.total_power = 1;
%         P.skip.eccentricity = 1;
%         P.skip.mean_angle = 1;
%         P.visible = 0;
%         P.HR = 0.98; % restored from before
%         PM = ICOS_Model6(P);
        %%
        res_a = res(i);
        res_a.Lenses = P.Lenses;
        res_a.Lens_Space = P.Lens_Space;
        res_a.detector_spacing = ds;
        res_a.theta = theta2;
        res_a.r_d = s*r / tand(theta2);
        res_a.Ltot = res_a.ORL + res_a.L + sum(res_a.Lens_Space) ...
          + res_a.detector_spacing;
        % res_a.overlap = PM.Results.overlap;
        results = results+1;
        res2(results) = res_a;
        fprintf(1,'%d: (%d,%d) %s%s @ %.2f (%.1f deg)\n', results, i, j, Lens1, LTS{j}, x, theta);
        %%
%         P.visible = 1;
%         PM = ICOS_Model6(P);
%         title(sprintf('Result %d: (%d,%d)', results, i, j));
%         drawnow;
      else
        warning('MATLAB:HUARP:BadAngle', 'Rejecting solution for bad angle');
      end
    end
  end
end