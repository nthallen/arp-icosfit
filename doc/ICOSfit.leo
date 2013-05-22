<?xml version="1.0" encoding="utf-8"?>
<?xml-stylesheet ekr_test?>
<leo_file>
<leo_header file_format="2" tnodes="0" max_tnode_index="0" clone_windows="0"/>
<globals body_outline_ratio="0.5" body_secondary_ratio="0.5">
	<global_window_position top="50" left="50" height="500" width="700"/>
	<global_log_window_position top="0" left="0" height="0" width="0"/>
</globals>
<preferences/>
<find_panel_settings/>
<vnodes>
<v t="ntallen.20130522094826.2380"><vh>@chapters</vh></v>
<v t="ntallen.20130522094826.2379" a="E"><vh>ICOSfit.exe</vh>
<v t="ntallen.20130522100930.1248" a="E"><vh>ToDo</vh>
<v t="ntallen.20130522124837.1267" a="E"><vh>Cleanup</vh>
<v t="ntallen.20130522124837.1268"><vh>Deprecated Configurations</vh></v>
<v t="ntallen.20130522124837.1269"><vh>ICOSsum.dat</vh></v>
<v t="ntallen.20130522124837.1270"><vh>ICOSconfig.m</vh></v>
</v>
<v t="ntallen.20130522100930.1249"><vh>Realtime Support</vh></v>
<v t="ntallen.20130522100930.1250"><vh>GSL Integration</vh></v>
<v t="ntallen.20130522100930.1251"><vh>Linking parameters</vh></v>
<v t="ntallen.20130522100930.1252"><vh>Parallel Operations</vh></v>
<v t="ntallen.20130522100930.1255"><vh>Documentation</vh></v>
</v>
<v t="ntallen.20130522100930.1253"><vh>Strategy</vh></v>
<v t="ntallen.20130522100930.1254"><vh>Troubleshooting</vh></v>
</v>
<v t="ntallen.20130522100930.1245"><vh>Matlab</vh>
<v t="ntallen.20130522100930.1246"><vh>etln_fit</vh>
<v t="ntallen.20130522100930.1247"><vh>etln_fit2</vh></v>
</v>
</v>
</vnodes>
<tnodes>
<t tx="ntallen.20130522094826.2379">@nocolor</t>
<t tx="ntallen.20130522094826.2380"></t>
<t tx="ntallen.20130522100930.1245">@nocolor</t>
<t tx="ntallen.20130522100930.1246">etln_fit:
    From processing CO2 data from Carbon isotopes
    Tau1 ~= 0.1
    Tau2 ~= 0.3
    [ ] I would like to reject if Tau1 drops below a minimal threshold

    When interactively adjusting Tau1 and Tau2, my current heuristic is:
      0: Adjust Tau1 (smaller of the two) to be above zero
      1: Adjust Tau1 until first few points are somewhat flat
      2: Adjust Tau2 until last few points are somewhat flat
      3: Iterate from 1 until both conditions are met
    [ ] Is it possible to automate this?
    [ ] Could this be implemented as a non-linear fit on just those two params?
    [ ] How close to the optimal do we need to get for the full fit to converge?
    [ ] How close to the optimal do we need to get for the two-param fit to converge?
        Monte-carlo.
</t>
<t tx="ntallen.20130522100930.1247">Currently an investigation. However, current systems have lasers that are stable enough that this work may not be necessary.

etln_fit2:
    For some range of scans
      Load scan
      autofindpeaks
      monte carlo {
        pick arbitrary values for tau1,2 (X([5 7]))  (.1, .3)
        use linear fit for X([1:4 6])
          plot residuals
        use polynomial fit for Y(8:11)
        set Y(12) = 1
          plot first guess
        do full fit and record residual
        plot starting tau1, tau2, ending residual
      }
      </t>
<t tx="ntallen.20130522100930.1248"></t>
<t tx="ntallen.20130522100930.1249"></t>
<t tx="ntallen.20130522100930.1250"></t>
<t tx="ntallen.20130522100930.1251"></t>
<t tx="ntallen.20130522100930.1252"></t>
<t tx="ntallen.20130522100930.1253">ICOSfit is somewhat more and somewhat less than a non-linear least squares fit.
It includes a non-linear least squares fitting engine, but is limited to
hard coded functional forms. It also includes a framework of other inputs
and parameters specifically designed to efficiently process multiple scans.
For better modularity in realtime and parallel implementations, we need to
make sure we record the necessary inputs.

Inputs:
    Configuration file
    Baseline configuration file
    Etalon fit inputs for tuning rate
    Initial line positions from matchline (or previous fit)
    Fitting range determined from line positions
    Dynamically fixed/floated parameters and thresholds

Voigt lines have 4 true fit parameters:
    N: Number density or Concentration, Molecules*cm-3
    fine_position: minor position adjustment, cm-1
    doppler width: cm-1
    lorentzian width: cm-1

Unless the line is outside the signal region, N remains free. There are independent
configuration controls for the other three parameters to indicate whether or not
they can float, but in order for these to float, the line strength must also
be above the specified threshold. This addresses one of the common failure modes
when fitting weak lines--that the fit will wander off the line or the line
width will become ridiculously large.

So there are two additional parameters per line:
    S_thresh: Threshold
    fixed: whether the line is currently fixed
There is hysteresis in the evaluation of the line strength. If the strength falls below
S_thresh, the line is fixed. It must then rise above 4*S_thresh in order to be floated
again. Hence knowing the line strength and the threshold is not enough to determine
whether or not the line is currently fixed.

The threshold value can by dynamically adjusted if the line is deemed to have wandered
close to the edge of the sample region. This can occur at high pressure where a scaling
factor is used for line margin. This means that subsequent fits may be fixed when they
otherwise would not be. [Perhaps this test should be limited to cases where the line is
close to only one edge]

Resetting the threshold and restarting can be a successful strategy when the fit fails
</t>
<t tx="ntallen.20130522100930.1254">Just try restarting where the fit failed. This resets thresholds to configuration values.</t>
<t tx="ntallen.20130522100930.1255">ICOSsum.dat format
</t>
<t tx="ntallen.20130522124837.1267"></t>
<t tx="ntallen.20130522124837.1268">Remove unnecessary cruft:
    Configuration Options:
        PTFile: Obsolete input file format
        PandTFile: Obsolete input file format
        Binary: Make this the default
        ASCII: Drop this format
        MinimumFringeSpacing
        DSFRLimits
        TuningRate

</t>
<t tx="ntallen.20130522124837.1269">The output currently includes the following columns:
    1: Time (but always zero!)
    2: P
    3: T
    4: cal_flow (but not included in PTEFile, so always zero)
    5: inlet_flow (ditto)
    6: fileno
    7: RateS(?)
    8: normalized chi-squared
    9: dFN (parameter from an obsolete etalon fit)

Rewrite to:
    1: fileno
    2: P
    3: T
    4: normalized chi-squared</t>
<t tx="ntallen.20130522124837.1270">Institute a new variable:
    ICOSfit_format_ver = 2;</t>
</tnodes>
</leo_file>
