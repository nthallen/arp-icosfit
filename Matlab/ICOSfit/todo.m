% todo.m
% TODO {
%   etln_fit8: the GUI
% }
% 
% CR Data analysis procedures
% 
% Engineering Data {
%   Requires copying CRENG*.MAT into C:\Data\CR\<date> directory
%   and renaming to creng*.mat due to peculiarities in ne_setup's
%   handling of the time variable.
%   Then run ui_cr
% }
% 
% Scan Data {
%   Determine Mirror Loss {
%     Use ringview to get tau
%     Mirror Loss is l/(c*tau)
%   }
%   Baseline Determination {
%     Use loadscans to extract baselines from appropriate scans
%     Use svds to generate modes
%     use base5 to determine how many modes are required
%     use writeskewbase to output the result {
%       Need mirror loss, cavity length numbers, sample rate, tolerance
%       sample rate can now be obtained via waves_used if you specify
%       the ofile as a number instead of a string.
%     }
%   }
%   fitline generation {
%     use matchline5
%     requires fitline.dat {
%       9 columns should do it (maybe need a 10th 0 column for some older scripts?)
%       Use columns from HITRAN specified in ICOSfit web doc
%       (1 1 2 3 5 7 8 9 14 - except 14 is not always 14, but it's %03d format).
%     }
%     Should be able to get x range from qclicomp eventually!
%     Should be able to fill in
%       SampleRate
%       QCLI_Wave
%       Sigma
%       ICOSdir
%       PTFile
%     line positions need to compensate for skew
%   }
%   Fit etalon {
%     etln_fit7
%   }
%   icosfit {
%   }
%   analysis {
%   }
% }
% 
% Tuning Rate {
%   I0 = mean of offline region
%   I=Raw etalon-I0
%   Im = median filter
%   sum contiguous regions above median
%   discard regions that are significantly smaller than the
%   previous region (though we should probably complain if we go
%   too long without accepting a region)
%   When a region is accepted, locate the max value and assign
%   Pmax there.
%   
%   Interpolate Pmax, then anneal to make it smooth
%   Pmax/I is of form 1+F*sin(pi*f)^2 where F and f are smoothly
%   varying functions. My parabolic fit should work well on this
%   for finding peaks, which should initialize the f function.
%   We can pick off peak values here to initialize the F function
%   as well.
%   Interpolate and anneal.
% }
% 
% 030203.4 {
%   Summary {
%     10-minute baseline scan plus some hysteresis work with
% 	varying amounts of isotopic constituents.
%   }
%   Investigations {
%     Explore use of SVD in characterizing baseline evolution
% 	Look at evolution of harmonic components
% 	See if we can successfully identify cavity modes
%   }
%   4605 - 5000 Background spectrum from QCL on with purge
%   5001 - 5272 Turn off purge and add room air
%   base.mat includes scans over the entire range.
% 
%   Given U,S,V, I want to look at the individual modes and
%   interpret the results. Also, do a DFT on the modes and
%   see if I can pick out the cavity length.
% 
%   [ icos, etln ] = loadscans( [], [4605:5000] );
% 	Try:
%   [ icos, etln ] = loadscans( [], [4605:5000], [600:1750] );
%   [ icos, etln ] = loadscans( [], [4605:5000], [1:2000] );
% 
% 	generate a wavenumber scale:
% 
%   etlna = - e_anneal(etln) * .019805;
% 	Now for any scan i, we could plot something akin to
% 	amplitude vs wavenumber as:
%   plot( etlna(:,i), icos(:,i) );
% 
%   Nyquist criteria dictates that the longest feature we
%   can resolve is -.5 ./ diff(etlna(:,1)) cm
% 
%   The scale for the shortest feature is around
%   1/(max(etlna(:,1))-min(etlna(:,1)));
%   That's one cycle in the entire scan. While a standard
%   FFT doesn't give more resolution beyond that, it isn't
%   clear to me that you cannot get significantly more
%   resolution with a DFT.
% 
%   To look at particular frequency components
% 
%   f = [65:.1:80];
%   F = dft(etlna(:,1),icos(:,1),f);
%   plot(f,abs(F));
% 
%   Looks like there's a peak at 76.8 (scan 500) or 76.86 (scan 1).
%   Might be interesting to look at the temporal evolution of the
%   dft in that range during warmup. Also look at all the harmonics
%   up to the nyquist.
% 
%   Look at the percent magnitude of error with a quadratic fit
%   and with the SVD fit with n columns.
% 
%   base8(icos);
% }
