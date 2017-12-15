# arp-icosfit
## ICOSfit: Data analysis tools for ICOS instruments

ICOSfit is a non-linear least-squares fit algorithm originally developed for processing data from the [Anderson Group](https://www.arp.harvard.edu/)'s ICOS instrument. It has since been extended to support other optical configurations, such as Herriott Cells.

This archive consists of the following:
- **icosfit:** The core ICOSfit fitting program
- **Matlab/ICOSfit:** a suite of Matlab support utilities for working with ICOS data and ICOSfit
- **Matlab/queryHITRAN:** a tool for working with the HITRAN database
- **Matlab/RayTrace:** a library for simple 3-D ray tracing
- **Matlab/ICOS_Models:** a set of models developed using Matlab/RayTrace

## ICOSfit
ICOSfit is a C++ program for fitting spectra from scanning laser spectrometers. It specifically understands and compensates for the lineshape perturbations caused by integrated cavity output spectroscopy (ICOS). The [ICOSfit Manual](https://www.arp.harvard.edu/eng/das/manuals/icosfit.html) describes the data formats and the syntax required by the icosfit program. Most users will take advantage of the Matlab support utilites in the Matlab/ICOSfit directory. These can go a long way toward performing many of the routine operations in setting up your data analysis. That said, these programs have a large number of options, tweaks, knobs and buttons, and not all of them have been documented.
