
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

double getDeltaindMOB0CSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindMOBLCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindMOBQCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);

int CalculateindMOB0CSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * SumLeft, double * SumRight, double * ErrorLeft, double * ErrorRight);
int CalculateindMOBLCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * SumLeft, double * SumRight, double * ErrorLeft, double * ErrorRight);
int CalculateindMOBQCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * SumLeft, double * SumRight, double * ErrorLeft, double * ErrorRight);

int LogSplitInFileindMOBLCSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double SumLeft, double SumRight, double ErrorLeft, double ErrorRight, double errorTotal);
int LogSplitHeaderindMOBLCSplit(char * FileName);