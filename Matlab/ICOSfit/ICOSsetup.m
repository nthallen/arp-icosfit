% ICOSsetup
% Assumes base variable is defined or defines it to be 'ICOSout'
% Loads ICOSout.sum into fitdata and defines
% C               ICOS_debug      chi2            n_input_params  
% Chi             Nfit            col             n_line_params   
% En              P               cpci14          n_lines         
% Gair            S0 Scorr        delta dFN dnu   nu nu_P
% Ged             STdep           fitdata         nu_F0           
% Gedcalc         T               iso             nu_text         
% Gl              T0              lines           row             
% Glcalc          base            molwts          v               
% Gv              binary          n_base_params   
% Gvcalc          by_molecule     n_cols
%%%%%
% v identifies the columns of fitdata corresponding to each line
if ~exist( 'base', 'var') | length(base) == 0
  base = 'ICOSout';
end

ICOS_debug = 0;

save_base = base;
if exist( [ base '/ICOSsum.dat' ] )
  fitdata = load( [ base '/ICOSsum.dat' ] );
elseif exist( [ base '/ICOSsum.out' ] )
  fitdata = load( [ base '/ICOSsum.out' ] );
  ICOS_debug = 1;
  base = [ base '/../../..' ];
else
  error(['Cannot locate ICOSsum.dat or ICOSsum.out in ' base ]);
end
nu0 = 0;
if exist( [ base '/ICOSconfig.m' ], 'file' )
  addpath( base );
  eval( 'ICOSconfig;' );
  rmpath( base );
  if ~exist('lines','var')
    error('lines not defined in ICOSconfig.m');
  end
  n_line_params = lines(:,9);
  if any(n_line_params ~= 4)
    error('Cannot handle different line types!');
  end
  n_line_params = 4;
  n_lines = size(lines,1);
  if ~exist('CavLen','var')
    CavLen = 76.8;
  end
  if nu0 == 0
    if exist('n_abs_params', 'var')
      if size(lines,1) > 0
        nu0 = floor(lines(1,3));
      else
        nu0 = 0;
      end
    else
      nu0 = 0;
      n_abs_params = 0;
    end
  elseif ~exist('n_abs_params', 'var')
    n_abs_params = 0;
  end
    
  if ~exist('n_abs_line_params', 'var')
    n_abs_line_params = 0;
  end
else
  linefile = [ base '/fitline.dat' ];
  if ~exist( linefile, 'file')
    linefile = 'fitline.dat';
    if ~exist(linefile, 'file')
      linefile = '/usr/CR/fitline.dat';
      if ~exist(linefile,'file')
        error 'Cannot locate fitline.dat'
      end
    end
  end
  lines = load(linefile);
  n_input_params = 8;
  n_base_params = 5;
  n_lines = size(lines,1);
  n_abs_params = 0;
  n_abs_line_params = 0;
end

n_cols = n_input_params + 2 * ...
  ( n_base_params + n_abs_params + ...
    n_lines * (n_line_params + n_abs_line_params) );
if n_cols ~= size(fitdata,2)
  disp(n_cols);
  disp(size(fitdata,2));
  error('ICOSconfig values don''t agree with fitdata columns');
end

if n_input_params == 4
  cpci14 = fitdata(:,1);
  chi2 = fitdata(:,4);
else
  cpci14 = fitdata(:,6);
  chi2 = fitdata(:,8);
end
if n_abs_params
  nu_F0 = fitdata(:,n_input_params+n_base_params+1);
  dFN = fitdata(:,9);
else
  nu_F0 = fitdata(:,9);
  dFN = fitdata(:,10);
end
base = save_base;
clear save_base;

v = ([1:n_lines]-1)*(n_line_params+n_abs_line_params) + ...
  n_input_params + n_base_params + n_abs_params + n_abs_line_params + 1;
row = ones(1,length(v));
col = ones( size(fitdata,1), 1 );
Nfit = fitdata(:,v+2);
P_vec = fitdata(:,2);
P = P_vec * row; % P in torr
T_vec = fitdata(:,3);
T = T_vec * row;
% C = 2.68675e19 * (P / 760.) ./ (T / 273.15); % [M]
C = 2.685e19 * (P / 760.) ./ (T / 273.15); % [M] Sayres' number
C(C==0) = NaN;
Chi = Nfit ./ C;

% Line definitions
S0 = lines(:,4)'; % line strength cm-1/(mol*cm-2)
En = lines(:,6)'; % Lower state energy in cm-1
Gair = lines(:,5)'; % air-broadend halfwidth at Tref, Pref
nu = lines(:,3)';
delta = lines(:,8)';
nu_P = col*nu + col*delta.*P/760.; % pressure-corrected line positions
iso = (lines(:,1)*10 + lines(:,2))';
nu_text = {};
if exist('by_molecule','var') & by_molecule > 0
  names = isovals(iso,'name');
  for i=1:n_lines
    nu_text{i} = sprintf('%d: %.4f %s',i, nu(i),names{i});
  end
  nu_text = char(nu_text);
else
  by_molecule = 0;
  % nu_text = num2str( nu', '%.4f');
  for i=1:n_lines
    nu_text{i} = sprintf('%d: %.4f',i, nu(i));
  end
  nu_text = char(nu_text);
end
STdep = lines(:,7)';
T0 = 296.;
Cfact = 1.438789;
Boltzfact = exp(Cfact*(col*En).*(T-T0)./(T0*T));
StimEmis = (1-exp(-Cfact*(col*nu)./T))./(1-exp(-Cfact*(col*nu)/T0));
% StimEmis = 1;
Qfact = (T0./T).^1.5;
Scorr = (col*S0) .* Qfact .* Boltzfact .* StimEmis; % correct line strength
if length(iso) == 0
  molwts = zeros(size(iso));
else
  molwts = isovals(iso,'weight')';
end
Ged = fitdata(:,v+1);
Gl = fitdata(:,v+3);
dnu = fitdata(:,v);
% Gedcalc = 3.581e-7 * (col*nu) .* sqrt(T./(col*molwts*log(2)));
Gedcalc = 4.30213e-7 * (col*nu) .* sqrt(T./(col*molwts));
Glcalc = ((T0./T).^(col*STdep)) .* (col*Gair) .* (P/760.);
Gv = 0.5346*Gl + sqrt( 0.2166 * Gl.*Gl + .69315*Ged.*Ged);
Gvcalc = 0.5346*Glcalc + sqrt( 0.2166 * Glcalc.*Glcalc + .69315*Gedcalc.*Gedcalc);