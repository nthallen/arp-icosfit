%%
load('exexparam3save.mat');
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
  P.R1 = res(i).R1;
  P.R2 = res(i).R2;
  P.r1 = 3*2.54/2;
  P.r2 = res(i).D2*2.54/2;
  P.mirror_spacing = res(i).L;
  P.y0 = res(i).Rr2;
  P.dy = -res(i).ORd2;
  P.dz = res(i).ORs2;
  % P.CT2 = res(i).CT2;
  P.HRC = res(i).RR1;
  P.Hr = res(i).RD1*2.54/2;
  P.HCT = 0.4;
  P.herriott_spacing = res(i).ORL;

  P.stop_ICOS = 0;
  P.visible = 1;
  P.visibility = [0 0 0];
  P.focus = 1;
  P.ICOS_passes_per_injection = 100;
  P.max_rays = 3000;
  P.injection_scale = 1;
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
    P.Lenses = {'Lens1',''};
    P.Lens_Space = [0.2, 0];
    Li = 2;
    %% recalculate d and continue
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
    [x,theta,ds] = pick_lens_x(r,d,s,f,15,Lens.r);
    if isempty(x)
%       P.Lens_Space(Li) = 1;
%       PM = ICOS_Model6(P);
%       title(sprintf('%d: Failure', i));
%       pause;
    else
      P.Lens_Space(Li) = x;
      P.detector_spacing = ds;
      res_a = res(i);
      res_a.Lenses = P.Lenses;
      res_a.Lens_Space = P.Lens_Space;
      res_a.detector_spacing = ds;
      res_a.theta = theta;
      res_a.r_d = s*r / tand(theta);
      res_a.Ltot = res_a.ORL + res_a.L + sum(res_a.Lens_Space) ...
        + res_a.detector_spacing;
      results = results+1;
      res2(results) = res_a;
      fprintf(1,'%d: %s @ %.2f (%.1f deg)\n', i, LTS{j}, x, theta);
      %PM = ICOS_Model6(P);
      %title(sprintf('%d: Success', i));
      %drawnow; shg;
    end
  end
end