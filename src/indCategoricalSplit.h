
#ifndef LEFT
#define LEFT      0x01
#endif

#ifndef RIGHT
#define RIGHT     0x00
#endif

#ifndef TREATMENT
#define TREATMENT    1
#endif

#ifndef CONTROL
#define CONTROL      0
#endif

#ifndef EVENT
#define EVENT        1
#endif

#ifndef NOEVENT
#define NOEVENT      0
#endif

#ifndef uint
typedef unsigned int  uint;
#endif

double getDeltailCategoricalSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userY, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs, int nbrClass, int nbrObsMinInControlClass, int nbrObsMinInTreatmentClass);
int CalculateCategoricalSplit(double * dataY, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int nbrClass, int nbrObsMinInControlClass, int nbrObsMinInTreatmentClass, int * nbrLeft, int * nbrRight, double * TauMaxLeft, double * TauMaxRight,int * maxClassLeft, int * maxClassRigth);
int CalculateCategoricalOneClassSplit(double * dataY, double * dataG, char * indCut, int NodeNbrObrservation, char LeftOrRight, int TheClassI, int nbrObsMinInEachClass, int * nbrClass, double * Mean);
int LogSplitInFileCategoricalSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double TauMaxLeft, double TauMaxRight, double errorTotal,int maxClassLeft, int maxClassRigth);
int LogSplitHeaderCategoricalCSplit(char * FileName);

