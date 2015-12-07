classdef ICOS_sr_search < handle
  % ICOS_sr_search
  % Class to explore alignments based on the sr optimization.
  % See Also:
  %   ICOS_sr_search.ICOS_sr_search
  %   ICOS_sr_search.enumerate
  %   ICOS_sr_search.design_tolerance
  %   ICOS_sr_search.build_tolerance
  %   ICOS_sr_search.draw_summary
  %   ICOS_sr_search.explore_build_tolerance
  %   ICOS_sr_search.validate_interleave_tolerance
  %   ICOS_sr_search.validate_build_tolerance
  %   ICOS_sr_search.focus
  %   ICOS_sr_search.savefile
  
  properties (SetAccess = immutable)
    SRopt
  end
  
  properties (Hidden = true)
    m_min % This spot can overlap
    m_max % Can't handle more spots than this
    SolInc % Just a preallocation parameter
    RawSols % The number of elements in Summary
    NSol % The solution number we're currently working on
    rdtanth % detector radius times tan of acceptance angle
    Rmax % how far out to plot
  end
  
  properties (SetAccess = private)
    Summary % Where we store the possible solutions
  end
  
  methods
    function SR = ICOS_sr_search(varargin)
      % SR = ICOS_sr_search(options)
      % options are keyword, value pairs. Keywords may include
      % any field of the SRopt property, including:
      %   C: Laser coherence length or equivalent, cm
      %   B: Laser beam diameter, cm
      %   rd: Target radius at detector, cm
      %   th: Target angle at detector, degrees
      %   L: ICOS cell length, cm
      %   Rtolerance: Tolerance factor on mirror radius
      %   r_limit: Maximum r1, cm (based on mirror diameter)
      %   Rw1: Target Rw1, cm
      SR.SRopt.C = 3000; % Laser coherence length or equivalent, cm
      SR.SRopt.B = .4; % Beam diameter, cm
      SR.SRopt.rd = .09; % Target radius at detector, cm
      SR.SRopt.th = 14.5; % Target angle at detector, degrees
      SR.SRopt.L = 50; % ICOS Cell length, cm
      SR.SRopt.Rtolerance = 0.01; % Tolerance factor on radius of curvature
      SR.SRopt.r_limit = (2.5*2.54/2 - SR.SRopt.B); % Maximum r1, cm (based on mirror diameter)
      SR.SRopt.Rw1 = 1.1*SR.SRopt.B/2; % Target Rw1, cm
      SR.SRopt.mnc = 'sr';
      for i=1:2:length(varargin)-1
        fld = varargin{i};
        if isfield(SR.SRopt, fld)
          SR.SRopt.(fld) = varargin{i+1};
        % elseif isfield(IS.ISopt,fld)
        %   IS.ISopt.(fld) = varargin{i+1};
        else
          error('MATLAB:HUARP:badopt', 'Invalid option: "%s"', fld);
        end
      end
      SR.m_min = ceil(SR.SRopt.C/(2*SR.SRopt.L)); % This spot can overlap
      SR.m_max = floor(pi/asin(SR.SRopt.B/(2*SR.SRopt.r_limit)));
      SR.SolInc = 40; % Just a preallocation parameter
      SR.RawSols = SR.SolInc;
      SR.NSol = 0;
      SR.rdtanth = SR.SRopt.rd*tand(SR.SRopt.th);
    end
    
    function enumerate(SR)
      % SR.enumerate
      % Identifies all interleave patterns that meet both the focus and
      % overlap criteria.
      Summary(SR.RawSols) = ...
        struct('RLmin',[],'r_max',[],'R1',[],'R2',[],'m',[],'k',[], ...
        'Phi',[], 'Phi_n', [], 'Phi_p', [], 'Phi_N', [], 'Phi_P', [], ...
        'Phi_OK', [], 'dL', [], 'dBL', [], 'sel', []);
      for m=SR.m_min:SR.m_max
        for k = 1:floor(m/2);
          dN = evaluate_interleave(m,k);
          if isempty(dN)
            continue;
          end
          phi = pi/m;
          Phi = k*phi;
          R1pole = SR.SRopt.L/sin(Phi)^2;
          if R1pole > 2000
            R1 = interp1([0 2000],[SR.SRopt.L 2000],1:2000)';
            R2 = SR.SRopt.L*(R1-SR.SRopt.L)./(sin(Phi)^2*R1-SR.SRopt.L);
            RR = []; % L/(1+cos(Phi));
            SR.Rmax = 2000;
          else
            R1a = interp1([0 1001],[SR.SRopt.L R1pole],1:1000)';
            R1b = interp1([0 1000],[R1pole 2000],1:1000)';
            R1 = [R1a;R1pole;R1b];
            R2 = [ ...
              SR.SRopt.L*(R1a-SR.SRopt.L)./(sin(Phi)^2*R1a-SR.SRopt.L); ...
              NaN; ...
              SR.SRopt.L*(R1b-SR.SRopt.L)./(sin(Phi)^2*R1b-SR.SRopt.L)];
            RR = SR.SRopt.L./(1+cos(Phi)*(-1));
            SR.Rmax = min(RR*1.1,2000);
            RR = RR(RR<2000);
          end
          R2(R2 > 2000 | R2 < -2000) = NaN;
          % fprintf(1,'Break at R1 = %.1f\n', L/sin(Phi)^2);
          r_min = SR.SRopt.B/(2*sin(pi/m));
          r_max = sqrt(SR.rdtanth* ...
            sqrt(((R2-SR.SRopt.L).*R1.^2*SR.SRopt.L)./ ...
              ((R1-SR.SRopt.L).*(R1+R2-SR.SRopt.L))));
          r_max(~isnan(r_max)) = min(r_max(~isnan(r_max)),SR.SRopt.r_limit);
          rRR = interp1(R1,r_max,RR,'linear','extrap');
          % This (commented) calculation is wrong, but I haven't figured out why yet
          %   k2 = ((R1-L).*(R1.*R2-L))./((R2-L).*R1.^2.*L);
          %   k2(k2<0) = NaN;
          %   s1 = r_max.*sqrt(k2);
          r2 = abs(r_max.*cos(Phi).*R2./(R2-SR.SRopt.L));
          w2 = abs(r_max*sin(Phi)*cos(Phi).*R2./(R2-SR.SRopt.L));
          s1 = w2/SR.SRopt.L;
          
          % Rw1 = 1.1*B/2;
          RL = abs(SR.SRopt.Rw1./s1);
          rRL = interp1(R1,RL,RR,'linear','extrap');
          vok = r_max > r_min & r_max >= r2;
          vnok = ~vok;
          RLok = RL;
          RLok(vnok) = NaN;
          RLnok = RL;
          RLnok(vok) = NaN;
          
          SR.NSol = SR.NSol+1;
          if SR.NSol > SR.RawSols
            SR.RawSols = SR.RawSols + SR.SolInc;
            Summary(SR.RawSols).RLmin = [];
          end
          
          Summary(SR.NSol).RLmin = nanmin(RLok);
          mi = find(RLok == Summary(SR.NSol).RLmin);
          phi_min = asin(SR.SRopt.B./(2*r_max(mi)));
          Phi_n = Phi + (phi_min-phi)/dN(2);
          Phi_p = Phi + (phi_min-phi)/dN(1);
          Summary(SR.NSol).R1 = R1(mi);
          Summary(SR.NSol).R2 = R2(mi);
          Summary(SR.NSol).r_max = r_max(mi);
          Summary(SR.NSol).m = m;
          Summary(SR.NSol).k = k;
          Summary(SR.NSol).Phi = Phi;
          Summary(SR.NSol).Phi_n = Phi_n;
          Summary(SR.NSol).Phi_p = Phi_p;
          Summary(SR.NSol).sel = 0;
          SP.R1 = R1(mi);
          SP.R2 = R2(mi);
          SP.L = SR.SRopt.L;
          SP.Rw1 = SR.SRopt.Rw1;
          SP.RL = SP.Rw1/s1(mi);
          Res = exparam(SP);
          check_params(SR.NSol, Res);
        end
      end
      %
      SR.Summary = Summary(1:SR.NSol);
    end
    
    function design_tolerance(SR)
      % SR.design_tolerance
      % Calculate range of length values L that could be required
      % to maintain overlap criteria over the range of possible
      % mirror radii of curvature allowed by Rtolerance.
      %
      % See Also:
      %   ICOS_sr_search.build_tolerance
      for i=1:length(SR.Summary)
        R1 = SR.Summary(i).R1 * (1+SR.SRopt.Rtolerance*[-1,1,-1,1]);
        R2 = SR.Summary(i).R2 * (1+SR.SRopt.Rtolerance*[-1,-1,1,1]);
        SR.Summary(i).Phi_OK = false;
        Phi_a = SR.SRopt.L.*(R1+R2-SR.SRopt.L)./(R1.*R2);
        if all(Phi_a >= 0 & Phi_a <= 1)
          Phi_b = asin(sqrt(Phi_a));
          Phi = minmax(Phi_b);
          SR.Summary(i).Phi_N = Phi(1);
          SR.Summary(i).Phi_P = Phi(2);
          if Phi(1)>= SR.Summary(i).Phi_n && Phi(2) <= SR.Summary(i).Phi_p
            SR.Summary(i).Phi_OK = true;
            SR.Summary(i).dL = [0;0];
          else
            v = Phi_b > SR.Summary(i).Phi_p;
            if any(v)
              % L necessary to reach Phi_p with R1(v), R2(v)
              det = (R1(v)+R2(v)).^2 - 4*R1(v).*R2(v).*sin(SR.Summary(i).Phi_p)^2;
              Lp = ([1;1]*(R1(v)+R2(v)) + [-1;1]*sqrt(det))/2;
              Lpr = Lp([1;1]*det > 0 & Lp > 0)';
            else
              Lpr = [];
            end
            v = Phi_b < SR.Summary(i).Phi_n;
            if any(v)
              det = (R1(v)+R2(v)).^2 - 4*R1(v).*R2(v).*sin(SR.Summary(i).Phi_n)^2;
              Ln = ([1;1]*(R1(v)+R2(v)) + [-1;1]*sqrt(det))/2;
              Lnr = Ln([1;1]*det > 0 & Ln > 0)';
            else
              Lnr = [];
            end
            SR.Summary(i).dL = minmax([SR.SRopt.L Lpr Lnr])' - SR.SRopt.L;
          end
        else
          SR.Summary(i).dL = [-2000;2000];
        end
      end
    end
    
    function draw_summary(SR)
      % SR.draw_summary;
      % Plot three axes showing RL_min, r1_max and R1+R2 vs
      % the advance angle, Phi.
      figure;
      ax = [ nsubplot(3,1,1), nsubplot(3,1,2), nsubplot(3,1,3)];
      Phi = [SR.Summary.Phi];
      RLmin = [SR.Summary.RLmin];
      Phi_OK = [SR.Summary.Phi_OK];
      r_max = [SR.Summary.r_max];
      plot(ax(1),Phi(~Phi_OK), RLmin(~Phi_OK), 'r.',Phi(Phi_OK), RLmin(Phi_OK), 'b*');
      plot(ax(2),Phi(~Phi_OK), r_max(~Phi_OK), 'r.',Phi(Phi_OK), r_max(Phi_OK), 'b*');
      plot(ax(3),[SR.Summary.Phi], [SR.Summary.R1], '*', [SR.Summary.Phi], [SR.Summary.R2], '+');
      title(ax(1),sprintf('L=%d', SR.SRopt.L));
      set(ax(1),'XTickLabel',[]);
      set(ax(2),'XTickLabel',[],'YAxisLocation','Right');
      ylabel(ax(1),'RL_{min} cm');
      ylabel(ax(2),'r_{1,max} cm');
      ylabel(ax(3),'R_{1,2}');
      grid(ax(3),'on');
      xlabel(ax(3),'Phi radians');
      linkaxes(ax, 'x');
    end
    
    function build_tolerance(SR)
      % SR.build_tolerance
      % Calculate range of length values over which L can vary
      % while maintaining overlap criteria.
      %
      % See Also:
      %   ICOS_sr_search.design_tolerance
      for i=1:length(SR.Summary)
        Phi = SR.Summary(i).Phi;
        R1 = SR.Summary(i).R1;
        R2 = SR.Summary(i).R2;
        det = (R1+R2).^2 - 4*R1.*R2.*sin([SR.Summary(i).Phi_n SR.Summary(i).Phi_p]).^2;
        if all(det > 0)
          Lp = ((R1+R2) + [-1;1]*sqrt(det))/2;
          Lpr = Lp(Lp>0);
          SR.Summary(i).dBL = minmax(Lpr')' - SR.SRopt.L;
        end
      end
    end
    
    function explore_build_tolerance(SR, plot_num)
      % SR.explore_build_tolerance([plot_num]);
      % plot_num 1 (default) is build tolerance vs RLmin
      % plot_num 2 is build tolerance vs design tolerance
      RLmin = [SR.Summary.RLmin];
      dL = [SR.Summary.dL];
      ddL = diff(dL);
      dBL = [SR.Summary.dBL];
      ddBL = diff(dBL);
      sel = [SR.Summary.sel];
      v = ddL < 2000;
      v1 = ddL < 2000 & sel;
      m = [SR.Summary.m];
      k = [SR.Summary.k];
      if nargin < 2 || plot_num == 1
        figure;
        h = plot(RLmin(v),ddBL(v),'.', RLmin(v1),ddBL(v1),'or');
        xlabel('RLmin cm');
        ylabel('Build \Delta{L} cm');
        title(sprintf('%s: Build Tolerance', SR.SRopt.mnc));
        hdt = datacursormode;
        set(hdt,'UpdateFcn', ...
          {@ICOS_sr_search.data_cursor_text_func,RLmin,ddBL,m,k,ddL,ddBL,RLmin});
        set(h,'buttondownfcn', @(s,e) ICOS_sr_search.ex_bdf(s,e,SR,RLmin,ddBL));
        set(gca, 'buttondownfcn', @(s,e) ICOS_sr_search.ex_bdf(s,e,SR,RLmin,ddBL));
      else
        figure;
        h = plot(ddL(v), ddBL(v),'.', ddL(v1), ddBL(v1), 'or');
        xlabel('Design \Delta{L} cm');
        ylabel('Build \Delta{L} cm');
        title(sprintf('%s: Build Tolerance', SR.SRopt.mnc));
        hdt = datacursormode;
        set(hdt,'UpdateFcn', ...
          {@ICOS_sr_search.data_cursor_text_func,ddL,ddBL,m,k,ddL,ddBL,RLmin});
        set(h,'buttondownfcn', @(s,e) ICOS_sr_search.ex_bdf(s,e,SR,ddL,ddBL));
        set(gca, 'buttondownfcn', @(s,e) ICOS_sr_search.ex_bdf(s,e,SR,ddL,ddBL));
      end
    end
    
    function validate_interleave_tolerance(SR,m,k)
      % SR.validate_interleave_tolerance(m,k)
      % Check that the Phi_n and Phi_p values actually define the limits.
      i = find([SR.Summary.m]==m & [SR.Summary.k]==k);
      if length(i) ~= 1
        return;
      end
      Phi = linspace(SR.Summary(i).Phi_n,SR.Summary(i).Phi_p,11);
      for j=1:length(Phi)
        draw_spots(SR.SRopt.r_limit+SR.SRopt.B, SR.Summary(i).r_max, ...
          2*Phi(j), m, SR.SRopt.B);
        shg;
        pause;
      end
    end
    
    function validate_build_tolerance(SR,m,k)
      % SR.validate_build_tolerance(m,k)
      % The build tolerance calculation seeks to determine how much the
      % cell length can vary without running into overlap trouble. This
      % method seeks to validate that claim by calculating alignments for a
      % number of points over the range and displaying spot patterns.
      i = find([SR.Summary.m]==m & [SR.Summary.k]==k);
      if length(i) ~= 1
        return;
      end
      Llims = SR.SRopt.L + SR.Summary(i).dBL;
      % L = linspace(Llims(1), Llims(2), 21);
      L = [Llims(1) SR.SRopt.L Llims(2)];
      R1 = SR.Summary(i).R1;
      R2 = SR.Summary(i).R2;
      Phi_a = L.*(R1+R2-L)./(R1.*R2);
      if all(Phi_a >= 0 & Phi_a <= 1)
        Phi = asin(sqrt(Phi_a));
        % Phi = minmax(Phi_b);
        r_min = SR.SRopt.B/(2*sin(pi/m));
        r_max = sqrt(SR.rdtanth* ...
          sqrt(((R2-L).*R1.^2.*L) ./ ((R1-L).*(R1+R2-L))));
        r_max(~isnan(r_max)) = min(r_max(~isnan(r_max)),SR.SRopt.r_limit);
        % figure;
%         plot(L,r_max,L,ones(size(L))*r_min);
%         legend('r_{max}','r_{min}');
%         xlabel('L');
%         ylabel('r');
%         shg;
%         pause;
        for j=1:length(L)
          draw_spots(SR.SRopt.r_limit+SR.SRopt.B, SR.Summary(i).r_max, ...
            2*Phi(j), m, SR.SRopt.B);
          shg;
          pause;
        end
      else
        fprintf(1,'Some bad angles\n');
      end
    end
    
    function focus(SR, varargin)
      % SR.focus
      % For each selected configuration, generate an ICOS_search
      % model and run search_ICOS_RIM and search_focus2 methods
      SFopt.focus_visible = 1;
      for i=1:2:length(varargin)-1
        fld = varargin{i};
        if isfield(SFopt, fld)
          SFopt.(fld) = varargin{i+1};
        else
          error('MATLAB:HUARP:badopt', 'Invalid option: "%s"', fld);
        end
      end
      selected = find([SR.Summary.sel]);
      Sums = SR.Summary(selected);
      for i = 1:length(selected)
        ii = selected(i);
        mnc = sprintf('%s.%d_%d.%d',SR.SRopt.mnc,ii,Sums(i).m,Sums(i).k);
        IS_fname = sprintf('IS_%s.mat', mnc);
        if ~exist(IS_fname,'file')
          fprintf(1,'Focusing (%d,%d)\n', Sums(i).m, Sums(i).k);
          SP.R1 = Sums(i).R1;
          SP.R2 = Sums(i).R2;
          SP.L = SR.SRopt.L;
          SP.Rw1 = 1.1*SR.SRopt.B/2;
          phi = pi/Sums(i).m;
          Phi = Sums(i).k*phi;
          r_max = sqrt(SR.rdtanth*sqrt(((SP.R2-SP.L).*SP.R1.^2*SP.L)./((SP.R1-SP.L).*(SP.R1+SP.R2-SP.L))));
          if isnan(r_max)
            error('r_max %d is NaN', i);
          end
          r_max = min(SR.SRopt.r_limit,r_max);
          r2 = abs(r_max.*cos(Phi).*SP.R2./(SP.R2-SP.L));
          w2 = abs(r_max*sin(Phi)*cos(Phi).*SP.R2./(SP.R2-SP.L));
          s1 = w2/SP.L;
          SP.RL = SP.Rw1/s1;
          Res = exparam(SP);
          check_params(i, Res);
          IS = ICOS_search('mnc', mnc,'R1',SP.R1,'R2',SP.R2,'L',SP.L, ...
            'RR1',Res.RR1,'Rw1',SP.Rw1, ...
            'focus_visible', SFopt.focus_visible);
          %%
          IS.search_ICOS_RIM;
          %%
          IS.search_focus2('max_lenses',2);
        else
          fprintf('%s already exists, skipping\n', IS_fname);
        end
      end
    end
    
    function savefile(SR)
      % SR.savefile
      % Writes the SR object to a .mat file with the name derived
      % from the mnemonic in SR.SRopt.mnc. The filename is prefixed
      % with 'SR_'.
      % See also: ICOS_beam
      fname = sprintf('SR_%s.mat', SR.SRopt.mnc);
      save(fname, 'SR');
      fprintf(1, 'ICOS_sr_search saved to %s\n', fname);
    end
    
  end
  
  methods(Static)
    function output_txt = data_cursor_text_func(~,event_obj,...
        x,y,m,k,ddL,ddBL,RL)
      % Display the position of the data cursor
      % obj          Currently not used (empty)
      % event_obj    Handle to event object
      % output_txt   Data cursor text string (string or cell array of strings).
      
      pos = get(event_obj,'Position');
      % output_txt = {['X: ',num2str(pos(1),4)],...
      %    ['Y: ',num2str(pos(2),4)]};

      % If there is a Z-coordinate in the position, display it as well
      % if length(pos) > 2
      %   output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
      % end
      
      output_txt = {};
      i = find(pos(1)==x & pos(2)==y,1);
      if ~isempty(i)
        output_txt{end+1} = sprintf('(m,k) = (%d,%d)', m(i), k(i));
        output_txt{end+1} = sprintf('ddL = %.2f', ddL(i));
        output_txt{end+1} = sprintf('ddBL = %.2f', ddBL(i));
        output_txt{end+1} = sprintf('RL = %.2f', RL(i));
      end
    end
    
    function ex_bdf(src,~,SR,X,Y)
      srctype = get(src,'type');
      while ~strcmp(srctype,'axes')
        src = get(src,'parent');
        if isempty(src)
          fprintf(1,'No axes found\n');
          return;
        end
        srctype = get(src,'type');
      end
      set(gcf,'units','normalized');
      axesHandle  = get(src,'Parent');
      coordinates = get(axesHandle,'CurrentPoint');
      coordinates = coordinates(1,1:2);
      % fprintf(1,'Click at coordinates: %f, %f\n', coordinates(1), coordinates(2));
      rect = rbbox([coordinates 0 0]);
      % fprintf(1,'Final rect: [%.2f %.2f %.2f %.2f]\n', rect);
      xl = get(gca,'xlim');
      yl = get(gca,'ylim');
      ap = get(gca,'position');
      XL = interp1([ap(1) ap(1)+ap(3)],xl,[rect(1),rect(1)+rect(3)]);
      YL = interp1([ap(2) ap(2)+ap(4)],yl,[rect(2),rect(2)+rect(4)]);
      % fprintf(1,'Final rect: X: [%.2f, %.2f] Y: [%.2f, %.2f]\n', ...
      %   XL, YL);
      v = find(X >= XL(1) & X <= XL(2) & Y >= YL(1) & Y <= YL(2));
      for i=1:length(v)
        SR.Summary(v(i)).sel = ~SR.Summary(v(i)).sel;
      end
      sel = find([SR.Summary.sel]);
      h = get(gca,'children');
      if length(h) == 2
        m = get(h,'marker');
        hi = find(strcmp(m, 'o'));
        set(h(hi),'XData',X(sel),'YData',Y(sel));
      else
        hold on;
        plot(X(sel),Y(sel),'or');
        hold off;
      end
    end
  end
end
