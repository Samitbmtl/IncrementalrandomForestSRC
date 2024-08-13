
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

double getDeltaindQuadCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindQuad2CSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindMaxQCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);
double getDeltaindMaxMaxQCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs);

int CalculateindQuadCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta0Left, double * Beta0Right, double * Beta1Left, double * Beta1Right, double * Beta2Left, double * Beta2Right);

int LogSplitInFileindQuadCSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double Beta0Left, double Beta0Right, double Beta1Left, double Beta1Right, double Beta2Left, double Beta2Right, double errorTotal);
int LogSplitHeaderindQuadCSplit(char * FileName);