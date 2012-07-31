function f = rrplot( base, basefit, CPCI14 );
% rrplot( base, basefit, CPCI14 );
% base is the data directory
% basefit is the quadratic fit parameters
% CPCI14 is the file number to plot.
% returns the data from the fit file.
path = mlf_path( base, CPCI14 );
f = load(path);
% baseline = polyval(basefit,f(:,1));
% baseline = basefit(1)*f(:,1).*f(:,1) + basefit(2)*f(:,1) + basefit(3);
plot(f(:,1),f(:,[2:4]));
grid;
title(path);

