classdef opt_model_p
  % Parametrized optical model
  properties
    M % opt_model
    Results % Results array
  end
  
  properties (Access = protected)
    p1
    vals1
    p1i
    p2
    vals2
    p2i
    writerObj
  end
  
  methods
    function PM = opt_model_p(P,p1,vals1,p2,vals2,meshz)
      % PM = opt_model_p(P, I);
      % P: Struct of model parameters
      Res = PM.results_struct;
      if nargin < 2
        p1 = '';
        vals1 = 1;
      else
        Res.(p1) = [];
      end
      PM.p1 = p1;
      PM.vals1 = vals1;
      if nargin < 4
        p2 = '';
        vals2 = 1;
      else
        Res.(p2) = [];
      end
      PM.p2 = p2;
      PM.vals2 = vals2;
      summary(length(vals2),length(vals1)) = Res; %preallocate
      PM.writerObj = [];
      
      for i = 1:length(PM.vals1)
        PM.p1i = i;
        if ~isempty(PM.p1)
          P.(p1) = vals1(PM.p1i);
        end
        for j = 1:length(PM.vals2)
          PM.p2i = j;
          if ~isempty(PM.p2)
            P.(PM.p2) = vals2(PM.p2i);
          end
          PM.M = PM.P_model(P);
          PM.M.propagate;
          Res = PM.evaluate_endpoints(P);
          Res.max_rays = P.max_rays;
          if ~isempty(p1)
            Res.(p1) = P.(p1);
          end
          if ~isempty(p2)
            Res.(p2) = P.(p2);
          end
          summary(PM.p2i,PM.p1i) = Res;
          if PM.visualize(P)
            PM.record_visualization(P);
          end
        end
      end
      if PM.M.visible && isfield(P, 'avifile')
        close(PM.writerObj);
        PM.writerObj = [];
      end
      PM.Results = PM.extract_results(summary);
      PM.Results.x = PM.p1;
      PM.Results.y = PM.p2;
      if nargin > 5
        PM.plot_results(meshz);
      end
    end
    
    function drawn = visualize(PM, P)
      % drawn = PM.visualize(P)
      % Returns true if something was drawn
      drawn = false;
      if PM.M.visible
        PM.draw_iteration_title(P);
        drawn = true;
      end
    end
   
    function ttl = draw_iteration_title(PM, P)
      if isfield(P,'title')
        p1ttl = '';
        if ~isempty(PM.p1)
          fmtfld = [ 'fmt_' PM.p1 ];
          if isfield(P, fmtfld)
            fmt = P.(fmtfld);
            p1ttl = sprintf([': ' fmt ], P.(PM.p1));
          else
            p1ttl = sprintf(': %d', PM.p1i);
          end
        end
        p2ttl = '';
        if ~isempty(PM.p2)
          fmtfld = [ 'fmt_' PM.p2 ];
          if isfield(P, fmtfld)
            fmt = P.(fmtfld);
            p2ttl = sprintf([', ' fmt ], P.(PM.p2));
          else
            p2ttl = sprintf(' x %d', PM.p2i);
          end
        end
        ttl = [ P.title p1ttl p2ttl ];
        title(ttl);
      end
    end
    
    function record_visualization(PM, P)
      if isfield(P, 'avifile')
        if isempty(PM.writerObj)
          PM.writerObj = VideoWriter(P.avifile);
          open(PM.writerObj);
        end
        frame = getframe(gcf);
        writeVideo(PM.writerObj, frame);
      elseif isfield(P, 'pause') && P.pause;
        drawnow;
        shg;
        pause;
      else
        drawnow; shg;
      end
    end
    
    function plot_results(PM, meshz)
      if isempty(PM.p1)
        fprintf(1,'p1 undefined in plot_results\n');
        return;
      end
      if isempty(PM.p2)
        clf;
        plot(PM.Results.(PM.p1),PM.Results.(meshz),'.');
        xlabel(strrep(PM.p1,'_','\_'));
        ylabel(strrep(meshz,'_','\_'));
      else
        try
          clf;
          mesh(PM.Results.(PM.p1),PM.Results.(PM.p2),PM.Results.(meshz));
          xlabel(strrep(PM.p1,'_','\_'));
          ylabel(strrep(PM.p2,'_','\_'));
          zlabel(strrep(meshz,'_','\_'));
        catch ME
          fprintf(1, 'mesh failed or invalid args: %s\n', ME.message);
        end
      end
    end
    
    function Res = evaluate_endpoints(PM,~)
      % Res = PM.evaluate_endpoints(P);
      % This is the default that does not know what the parameter
      % structure holds, so it does not use it and simply
      % evaluates the last optic in the chain.
      % This can be overridden to give paramteric control.
      opt_n = length(PM.M.Optic);
      Res = PM.M.evaluate_endpoints(opt_n);
    end
    
    function PM = clean_results(PM)
      % PM = PM.clean_results;
      % Using PM.Results.n_rays, marks all all results as NaN
      % where n_rays is not maximal. The assumption is that
      % rays have escaped.
      if isfield(PM.Results,'n_rays')
        snr = size(PM.Results.n_rays);
        mnr = max(max(PM.Results.n_rays));
        subpar = PM.Results.n_rays < mnr;
        flds = fieldnames(PM.Results);
        for i=1:length(flds)
          fld = flds{i};
          if all(size(PM.Results.(fld)) == snr)
            PM.Results.(fld)(subpar) = NaN;
          end
        end
      end
    end
  end
  
  methods (Abstract, Static)
    P = props
    M = P_model(P)
    % P = opt_model_p.props;
    % Returns struct P identifying appropriate parameters for
    % the model.
  end
    
  methods (Static)
    function Res = results_struct
      Res = opt_model.results_struct;
    end
    
    function R2 = extract_results(R1)
      % R2 = extract_results(R1);
      % R1: array of structs
      % R2: struct of arrays
      flds = fieldnames(R1);
      for i=1:length(flds)
        fld = flds{i};
        R2.(fld) = reshape([R1.(fld)],size(R1,1),size(R1,2));
      end
    end
  end
end
