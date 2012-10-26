function rrfit( base, range, plotcode, xlim_in, ylim_in );
% rrfit( base [, range [, plotcode]] );
% Display ICOS fits found in the specified directory.
% plotcode specifies what sort of plot to generate
% plotcode == 0 (default) is the current plot du jour
% plotcode == 1 plot fit
% plotcode == 2 plot strength instead of fit
% plotcode == 3 plot strength and fit
% plotcode == 5 plot fit, res, and base as %
% plotcode == 10 show per-line error in mixing ratio
% plotcode == 11 show per-line error in number density
% plotcode == 12 show per-line error in line width gamma_ed
% plotcode == 20 show error in wavenumber scale
% plotcode == 21 show error in baseline fit
if nargin < 5
  ylim_in = [];
  if nargin < 4
    xlim_in = [];
    if nargin < 3
      plotcode = 0;
      if nargin < 2
        range = [];
      end
    end
  end
end
lines=[];
by_molecule=0;
ICOSsetup;

plotans = 0;
plotfit = 1;
plotstrength = 0;
plotbase = 0;
if plotcode == 0
  plotcode = 1;
end
if plotcode >= 10 && plotcode < 30
  % load the answerkey
  ICOSfit_cfg = load_ICOSfit_cfg;
  run = getrun(1);
  ak = load([ ICOSfit_cfg.Matlab_Path filesep run '/answerkey.mat' ]);
  aChi = ak.ppm(scannum)*1e-6*row;
  aN = aChi .* C;
  PowerScale = ak.mirrorloss * 1e-6 / (2 - ak.mirrorloss*1e-6);
  if plotcode < 20; plotans = 1; end
end

if plotcode == 1
  plotfit = 1;
  plotstrength = 0;
elseif plotcode == 2
  plotfit = 0;
  plotstrength = 1;
elseif plotcode == 3
  plotfit = 1;
  plotstrength = 1;
elseif plotcode == 30
    plotfit = 1;
    plotstrength = 0;
    plotbase = 1;
elseif plotcode == 10
  % answer plot using Mixing Ratio Error
  errlabel = '%err MR';
  errval = 100*((Chi./aChi)-1);
elseif plotcode == 11
  % answer plot using Number Density
  errlabel = '%err ND';
  errval = 100*((Nfit./aN)-1);
elseif plotcode == 12
  % answer plot using Ged
  errlabel = '%err \gamma_{ED}';
  errval = 100*(Ged./Gedcalc - 1);
end
% nu = lines(:,3)';
% delta = lines(:,8)';

if length(range)==0
  rows = [1:size(fitdata,1)]';
else
  if size(range,2) > 1
    range = range';
  end
  rows = unique(interp1( scannum, [1:length(scannum)]', range, 'nearest' ));
  rows = rows(isfinite(rows));
  % rows = find(scannum >= min(range) & scannum <= max(range));
end

figno = figure;
for i=rows'
  lasterr('');
  % try
  if ICOS_debug
    path = sprintf( '%s/%04d.dat', base, scannum(i));
  else
    path = mlf_path( base, scannum(i));
  end
  if exist(path,'file')
    f = load(path);
    
    lpos = nu + delta*P_vec(i)/760. - fitdata(i,v);
    if nu0 ~= 0
      lpos = lpos + nu_F0(i);
    end
%   lpos = nu + delta*P(i)/760.;
    lgd = fitdata(i,v+1);
    lst = Scorr.*CavLen.*Nfit./(Ged*sqrt(pi));
    lgl = fitdata(i,v+3);
    lwid = (lgd+lgl); % fitdata(i,v+1)+fitdata(i,v+3);
    nux = f(:,2);
    if min(nux) < 2
      nux = nux+nu_F0(i)+nu0; % This may change
    end
    if plotcode >= 30
        basep = interp1(nux,f(:,5)./f(:,5)*100,lpos,'nearest');
        fitp =  interp1(nux,f(:,4)./f(:,5)*100,lpos,'nearest');
        meanp = (basep+fitp)/2;
        % threw in maxp to see the position of small lines
        maxp = ones(size(lpos))*100;
    else
        basep = interp1(nux,f(:,5),lpos,'nearest');
        fitp =  interp1(nux,f(:,4),lpos,'nearest');
        meanp = (basep+fitp)/2;
        % threw in maxp to see the position of small lines
        maxp = ones(size(lpos))*max(f(:,5));
    end
    X = [ lpos lpos-lwid; lpos lpos+lwid ];
    % Y = [ basep meanp; fitp meanp ];
    Y = [ maxp meanp; fitp meanp ];
    if plotcode < 20
      residual = f(:,3)-f(:,4);
      reslbl = 'Fit Res';
    elseif plotcode == 20
      residual = ak.v(f(:,1))' - nux;
      reslbl = 'nu Res';
    elseif plotcode == 21
      if strcmp(ak.runname,'quadtocubic')
        V0 = [0 -60 1.333e5 0];
        V1 = [.011 -83 1.43e5 0];
        r = (scannum(i)-1)/200;
        V = (1-r)*V0 + r*V1;
        abase = polyval(V,ak.x') * PowerScale;
      else
        abase = ak.Power' * PowerScale;
      end
      
      residual = abase(f(:,1)) - f(:,5);
      reslbl = 'base Res';
    elseif plotcode == 30
        residual = (f(:,3)-f(:,4))./f(:,5)*100;
        reslbl = 'Fit Res (%)';
    else
      error('Invalid plotcode number');
    end
    if plotfit
      clf;
      sdev = sqrt(mean((residual).^2));
      ttltext = [ 'Scan:' num2str(scannum(i)) ...
            ' \sigma = ' num2str(sdev)  ' ' getrun '/' path ];
      if plotans; nsubplots=4; else nsubplots=3; end
      if plotbase==1; nsubplots=5; end
      nsubplot( nsubplots, 1, [ nsubplots 2 ], 1 );
      if plotcode == 30
        plot(nux,f(:,3)./f(:,5)*100,nux,f(:,4)./f(:,5)*100,X,Y,'r')
        yl=ylim;
        ylim([yl(1),100+diff(yl)*.025]);
       % xlim([2309.5,2309.7])
      else
        plot(nux,f(:,[3:5]), X, Y, 'r');
        %xlim([2309.5,2309.7])
      end
      %legend('raw','fit','base');
      set( gca, 'XDir', 'reverse','YAxisLocation','Right' );
      grid;
      xr = [min(nux) max(nux)];
      if length(xlim_in)
        xlim(xlim_in);
      %xlim([2309.5,2309.7])
      else
        xlim( xr + [-1 1]*.05*diff(xr));
       % xlim([2309.5,2309.7])
      end
      xl = xlim;
      %xl=[2309.5,2309.7];
      if length(ylim_in)
        ylim(ylim_in);
      end

      nsubplot( nsubplots, 1, 1, 1 );
      if plotans
        bar(lpos, errval(i,:));
        ylabel(errlabel);
        set(gca,'XTickLabel',[],'XDir','reverse');
        xlim(xl);
        title(ttltext);
        nsubplot( nsubplots, 1, 2, 1 );
      end

      plot(nux,residual); grid;
      set(gca,'XTickLabel',[],'XDir','reverse');
      %set(gca,'YAxisLocation','right');
      xlim(xl);
      ylabel(reslbl);
      if ~plotans; title(ttltext); end 
      
      if plotbase
          nsubplot( nsubplots, 1, 2, 1 );
          plot(nux,f(:,5),'r'); grid
          set(gca,'XTickLabel',[],'XDir','reverse','YAxisLocation','Right');
          xlim(xl);
          ylabel('Baseline');
          nsubplot( nsubplots, 1, 3, 1 );
          plot(nux,detrend(f(:,5)),'r'); grid
          set(gca,'XTickLabel',[],'XDir','reverse');
          xlim(xl);
          ylabel('Detrended Baseline');
      end
      
      addzoom;
      drawnow; shg;
      pause;
    end
    if plotstrength
      clf;
      nsubplot( 2, 1, 1, 1 );
      plot(nux,f(:,[3:5]), X, Y, 'r');
      grid;
      set(gca,'XTickLabel',[],'XDir','reverse');
      title(path);
      xl = xlim;

      nsubplot(2, 1, 2, 1 );
      % strength = fitdata(i,v+2);
      semilogy(lpos, lst(i,:), '+');
      ylim([1e-7 1e-2]);
      set(gca,'XDir','reverse');
      grid;
      xlim(xl);
      addzoom;
      drawnow; shg;
      pause;
    end
    if any(findobj('type','figure')==figno)
      figure(figno);
    else
      break;
    end
  end
end
