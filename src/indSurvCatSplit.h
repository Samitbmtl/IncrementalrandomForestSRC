
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

double getDeltaSurvCatSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime,double * userEvent,double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incTailType, int incMinObs, int nbrClass, int nbrObsMinInControlClass, int nbrObsMinInTreatmentClass);
int CalculateSurvCatSplit(double * dataT,double * dataC, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int nbrClass, int nbrObsMinInControlClass, int nbrObsMinInTreatmentClass, int incTailType, int * nbrLeft, int * nbrRight, double * TauMaxLeft, double * TauMaxRight,int * maxClassLeft, int * maxClassRigth);
int CalculateSurvCatOneClassSplit(double * dataT,double * dataC, double * dataG, char * indCut, int NodeNbrObrservation, char LeftOrRight, int TheClassI, int nbrObsMinInEachClass, double ** T, double ** S, int * ClassNbr);
int LogSplitInFileSurvCatSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double TauMaxLeft, double TauMaxRight, double errorTotal,int maxClassLeft, int maxClassRigth);
int LogSplitHeaderSurvCatSplit(char * FileName);