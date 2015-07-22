function P_o = render_model(res)
% render_model(res);
%   Invokes ICOS_Model6
% P = render_model(res);
%   Just returns the parameter structures.
P = ICOS_Model6.props;
if nargout > 0
  P_o(length(res)) = P;
end
for i=1:length(res)
  P.R1 = res(i).R1;
  P.R2 = res(i).R2;
  P.r1 = 3*2.54/2;
  P.r2 = res(i).D2*2.54/2;
  P.mirror_spacing = res(i).L;
  P.y0 = res(i).Rr2;
  P.dy = -res(i).ORd2;
  P.dz = res(i).ORs2;
  P.HRC = res(i).RR1;
  P.Hr = res(i).RD1*2.54/2;
  % P.HCT = 0.4;
  P.herriott_spacing = res(i).ORL;
  if isfield(res(i),'Lenses')
    P.Lenses = res(i).Lenses;
    P.Lens_Space = res(i).Lens_Space;
    P.detector_spacing = res(i).detector_spacing;
    P.focus = 1;
  else
    P.focus = 0;
  end

  P.stop_ICOS = 0;
  P.visible = 1;
  P.visibility = [];
  P.ICOS_passes_per_injection = 100;
  P.max_rays = 3000;
  P.injection_scale = 1;
  
  if nargout > 0
    P_o(i) = P;
  else
    figure;
    PM = ICOS_Model6(P);
  end
end
