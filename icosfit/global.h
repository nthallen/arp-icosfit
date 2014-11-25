#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED
#include "funceval.h"

class GlobalData_t {
  public:
	unsigned int BackgroundRegion[2];
	unsigned int SignalRegion[2];
	unsigned int ScanNumRange[2];
	unsigned int RestartAt;
	int PreserveOutput;
	const char *FitFunction;
	float MirrorLoss;
	float EtalonFSR;
  float EtalonFeedback;
	// float MinimumFringeSpacing;
	float TolerableDrift;
  float LeftLineMargin;
  float LeftLineMarginMultiplier;
  float RightLineMargin;
  float RightLineMarginMultiplier;
  float LineMarginHysteresis;
	// float DSFRLimits[2];
	float CavityLength;
  float CavityFixedLength;
	float Sigma;
	// float TuningRate;
	float SampleRate;
	float SkewTolerance;
	float DefaultTemp;
	unsigned short QCLI_Wave;
	int binary;
	const char *ICOSdir;
	const char *PTFile;
	int PTformat;
  int PTE_nu_F0_col;
  int PTE_MirrorLoss_col;
  int PTE_PowerParams_col;
  int PTE_Feedback_col;
	const char *BaselineFile;
	const char *LineFile;
	const char *OutputDir;
	const char *OutputFile;
	const char *LogFile;
	const char *MFile;
  const char *QTdir;
	int Verbosity;
	int N_Passes;
  int BaselineInput; // non-zero if column 3 is  a baseline vector
	func_abs_p absorb;
	struct {
	  float nu_F0;
    float MirrorLoss;
	} input;
  float ConvergenceStep;
  int ConvergenceCount;
  int MaxIterations;

	GlobalData_t();
};
extern GlobalData_t GlobalData;
#define SetGlobal(x,y) GlobalData.x = y
#define SetGlobalOnce(x,y) (GlobalData.x ? 1 : ((GlobalData.x = y), 0))
#define SetGlobalPair(x,y,z) GlobalData.x[0] = y; GlobalData.x[1] = z
#define SetNu0(x) func_line::nu0 = x
#define Nu0IsSet  (func_line::nu0 != 0)
extern void ICOS_init();
extern void ICOS_main();

#endif
