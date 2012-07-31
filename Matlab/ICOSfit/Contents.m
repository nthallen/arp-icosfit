% ICOS Viewing and Fitting
% 
% Raw Data:
%   icosview.m - Display logged ICOS and Ringdown data
%   ringview.m - Display just Ringdown data and bin two ways
%   icosnoise.m - Calculate noise statistics
%   rawview.m - Review logged RAW ringdown data and compare fits
%     logchi.m - called via fmins() by rawview.m
%     lvoffset.m - Display results of rawview.m
%   qcliview.m - display engineering data relevant to QCLI
% 
% ICOSfit support:
%   matchline5.m - Create ICOSfit configuration file
%   diagnose.m - View ICOSfit Details
%   mixlines.m - View ICOSfit mixing ratios
%   dispfix.m - View ICOSfit fix/float status
%   rrfit.m - Display individual fit files
%   rrcompare.m - Plots to compare fit files from two bases
%   writebase.m - Write baseline files
%   writeskewbase.m - Write de-skewed baseline files
% 
% Utilities:
%   fitfringe.m - Locate Diagnostic Etalon fringes
%   cpciload.m - locate and load specified CPCI14 file
%   loadbin.m - simple binary loader for CPCI14 files
%   writebin.m - produce icos-format binary files
%   ICOSsetup.m - Common utility for locating directories, etc.
%   mlf_path.m - Create path from CPCI14 number
%   peakfind.m - Generic peak finder. Used in matchline5, icosnoise
%   isovals.m - HITRAN values
%   fixcpci14.m - adjust for apparent time delays in reporting CPCI14
%   humlik.m - Port of humlik.for
% 
% Studies:
%   humlicek.for
%   humlik.for
%   cavity.m - Fit cavity etalon to residuals
% 	reschi.m
%   ettest.m - Study leading toward fit_fringe.m
%   etalon.m - Study leading toward fit_fringe.m
%   matchline2.m - Series of matchlines
%   matchline3.m - 
%   matchline4.m - Never actually in production
%   rfit.m - simple version of rrfit.m
%   rrfitx.m - working backwards from fit with synthetic spectra
%   rrplot.m - obsolete version of rrfit.m
%   rrtune.m - studies of tuning rate issues
%   rrtune2.m - studies of tuning rate issues
%   tune[2-6].m - Tuning rate studies
%   tunetest.m - Tuning rate study
%   t_error.m - Explore effects of errors in T on fits
% 
% 
% loadscans.m
% base4.m
% basefind.m
% basetest.m
% bfilter.m
% evolve.m
% hypetalon.m
% losttrig.m
% testfit.m
% testfit2.m
% testfringe.m
