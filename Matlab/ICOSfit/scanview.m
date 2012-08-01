function scanview( base, cpci14 );
% scanview( base, cpci14 );
% Display raw scan with etalon, optionally limited to a range
% of CPCI14 values.
% Attempting modification to support ringdown data as well
if nargin < 1
  base = '';
  has_cpci = 0;
elseif nargin == 1 & isnumeric(base)
  cpci14 = base;
  base = '';
  has_cpci = 1;
elseif nargin > 1
  has_cpci = 1;
else
  has_cpci = 0;
end
PT = load_mat_files('PT');
WaveSpecs = load_waves;
v = find(diff(PT.CPCI14)>0)+1; % index of new CPCI14 values
v = [ v(1)-1; v ];
if ~has_cpci
  cpci14 = PT.CPCI14;
  cpci14 = [min(cpci14(cpci14 > max(0,min(cpci14)))):max(cpci14)];
else
  % constrain to within the range of defined values
  cpci14 = cpci14(cpci14 > 0 & cpci14 > min(PT.CPCI14) & ...
    cpci14 <= max(PT.CPCI14));
end

% Now locate each specified cpci14 as an index into v, the
% unique PT.CPCI14 entry indexes
idx = v(ceil(interp1( PT.CPCI14(v), [1:length(v)], cpci14 )));
% idx is now an array as long as cpci14
% PT.CPCI14(idx) should be equal to cpci14 except where skipping
% occurs, and then it should be greater than cpci14.
%-----------------------------
% roris = Waves(PT.QCLI_Wave(idx)+1,7);
if exist('WaveSpecs','var')
  roris = ~[ WaveSpecs(PT.QCLI_Wave(idx)+1).ISICOS ]';
else
  roris = zeros(size(idx));
end

base = find_scan_dir(base);
binary = 1;
if size(cpci14,1) > 1; cpci14 = cpci14'; end
figno = figure;
for i = 1:length(cpci14)
  path = mlf_path( base, cpci14(i), '.dat');
  data_ok = 0;
  if binary
    fe = loadbin( path );
    data_ok = (length(fe)>0);
  else
    if exist(path,'file')
      fe = load(path);
      data_ok = 1;
    end
  end
  if data_ok
    nsamples = size(fe,1);
    figure(figno);
    clf;
    if roris(i) == 0
      nsubplot(2,1,1,1);
      plot([1:nsamples],fe(:,1));
      set(gca,'xticklabel',[]);
      title(sprintf('ICOS Scan %d: %s', cpci14(i), getrun ));
      % ylim([0 50000]);
      nsubplot(2,1,2,1);
      plot([1:nsamples],fe(:,2));
      set(gca,'YAxisLocation','right');
      xlabel('Samples');
      drawnow; shg;
      pause;
    else
      if any( ~(isnan(fe(:,6))|isnan(fe(:,1))))
        tau = .1 * fe(:,6) ./ log(fe(:,1));
        if size(fe,2) > 7
          nsubplot(2,1,1);
        end
        plot(tau);
        ylabel('Tau \mu secs');
        title(sprintf('Ringdown Scan %d: %s', cpci14(i), getrun ));
        ylim([12 14.5]);
        if size(fe,2) > 7
          set(gca,'XTickLabel', []);
          set(gca,'YAxisLocation', 'right');
          nsubplot(2,1,2);
          plot(fe(:,8));
          ylabel('Etalon');
        end
        xlabel('Offset');
        drawnow; shg;
        pause;
      end
    end
    
    if all(findobj('type','figure') ~= figno)
      break;
    end
  else
    fprintf( 1, 'cpci %d: %s not found or unreadable\n', cpci14(i), path );
  end
end
