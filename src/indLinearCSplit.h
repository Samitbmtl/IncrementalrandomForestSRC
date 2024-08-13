
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

double getDeltaindLinearCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindTSTCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindNSRCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindTST2CSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindMaxLCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);

int CalculateindLinearCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right);
int CalculateindTSTCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right);
//int CalculateindNSRCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right);
int CalculateindTST2CSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right);

int LogSplitInFileindLinearCSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double Beta1Left, double Beta1Right, double errorTotal);
int LogSplitHeaderindLinearCSplit(char * FileName);