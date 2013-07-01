#include "global.h"

GlobalData_t GlobalData;

GlobalData_t::GlobalData_t() {
  BackgroundRegion[0] = 10;
  BackgroundRegion[1] = 200;
  SignalRegion[0] = 350;
  SignalRegion[1] = 1750;
  ScanNumRange[0] = ScanNumRange[1] = 0;
  RestartAt = 0;
  PreserveOutput = 0;
  FitFunction = "func_tau";
  MirrorLoss = 180.e-6;
  EtalonFSR = 0.019805;
  MinimumFringeSpacing = 12.;
  TolerableDrift = .01; // cm-1
  CavityLength = 70.; // cm
  CavityFixedLength = 0.; // cm
  LeftLineMargin = RightLineMargin = .05; // cm-1
  LeftLineMarginMultiplier = RightLineMarginMultiplier = 8;
  LineMarginHysteresis = 1e-6;
  DSFRLimits[0] = .95;
  DSFRLimits[1] = 1.21;
  Sigma = 1000.;
  TuningRate = 0.;
  QCLI_Wave = 0;
  SampleRate = 0.;
  SkewTolerance = 1e-5;
  DefaultTemp = 22.56;
  binary = 1;
  ICOSdir = "Scans";
  PTFile = "INPUTFILE.txt";
  PTformat = 0;
  BaselineFile = 0;
  BaselineInput = 0;
  LineFile = "fitline.dat";
  OutputDir = "ICOSout";
  OutputFile = "ICOSsum.dat";
  LogFile = "ICOSfit.log";
  MFile = "ICOSconfig.m";
  QTdir = "/usr/local/share/QT";
  Verbosity = 0;
  absorb = 0;
  N_Passes = 0;
  ConvergenceStep = 1e-3;
  ConvergenceCount = 4;
  MaxIterations = 500;
}
